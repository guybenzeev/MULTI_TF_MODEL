from __future__ import annotations

import csv
from pathlib import Path


class Multi_Trajectory:
    def __init__(self, csv_name: str) -> None:
        data_dir = Path(__file__).resolve().parent.parent / "data"
        self.file_path = data_dir / csv_name

        if not self.file_path.is_file():
            raise FileNotFoundError(
                f"CSV file not found: {self.file_path}"
            )

        self.file_object = self.file_path.open(
            newline="",
            encoding="utf-8",
        )

        rows = list(csv.reader(self.file_object))

        self.header = rows[0] if rows else []
        self._metadata = self._parse_metadata(self.header)
        self.lines = rows[1:]
        self._current_index = 0

    @staticmethod
    def _parse_metadata(header: list[str]) -> dict[str, str]:
        metadata: dict[str, str] = {}

        # Event columns:
        # time, index, edit, protein, side
        #
        # Everything after those columns is metadata.
        for value in header[5:]:
            if "=" not in value:
                continue

            key, parsed_value = value.split("=", 1)
            metadata[key.strip()] = parsed_value.strip()

        return metadata

    def get_line(self, index: int) -> list[str]:
        return self.lines[index]

    def init_iterator(self) -> None:
        self._current_index = 0

    def next_line(self) -> list[str]:
        if not self.has_next():
            raise StopIteration("No more lines in the CSV file.")

        line = self.lines[self._current_index]
        self._current_index += 1
        return line

    def has_next(self) -> bool:
        return self._current_index < len(self.lines)

    def get_length(self) -> int:
        return int(self._metadata["length"])

    def get_num_sides(self) -> int:
        return int(self._metadata["num_sides"])

    def get_num_proteins(self) -> int:
        return int(self._metadata["num_proteins"])

    def get_runtime(self) -> float:
        return float(self._metadata["runTime"])

    @staticmethod
    def get_state_id(
        binding_type: str,
        protein: int,
        side: int,
        num_sides: int,
        num_proteins: int,
    ) -> int:
        """
        Reproduce the C++ State::getStateID schema.

        FREE:
            0

        BOUND_NS:
            1 + (protein - 1) * num_sides + (side - 1)

        BOUND_S:
            1 + num_sides * num_proteins
              + (protein - 1) * num_sides
              + (side - 1)
        """
        if binding_type == "FREE":
            return 0

        if protein <= 0 or protein > num_proteins:
            raise ValueError(
                f"Invalid protein {protein}. "
                f"Expected a value from 1 to {num_proteins}."
            )

        if side <= 0 or side > num_sides:
            raise ValueError(
                f"Invalid side {side}. "
                f"Expected a value from 1 to {num_sides}."
            )

        protein_index = protein - 1
        side_index = side - 1
        per_binding_type = num_sides * num_proteins

        if binding_type == "BOUND_NS":
            return (
                1
                + protein_index * num_sides
                + side_index
            )

        if binding_type == "BOUND_S":
            return (
                1
                + per_binding_type
                + protein_index * num_sides
                + side_index
            )

        raise ValueError(
            f"Unknown binding type: {binding_type}"
        )

    def applyEdit(
        self,
        line: list[str],
        list_state: list[int],
    ) -> None:
        if len(line) < 5:
            raise ValueError(
                "Expected CSV row with five columns: "
                "time,index,edit,protein,side"
            )

        _, index_text, edit_type, protein_text, side_text = line[:5]

        index = int(index_text)
        protein = int(protein_text)
        side = int(side_text)

        num_sides = self.get_num_sides()
        num_proteins = self.get_num_proteins()

        if index < 0 or index >= len(list_state):
            raise IndexError(
                f"Trajectory index {index} is outside the state "
                f"list of length {len(list_state)}."
            )

        if edit_type in {
            "BIND_NS",
            "SWITCH_SIDE",
            "UNBIND_S",
        }:
            list_state[index] = self.get_state_id(
                binding_type="BOUND_NS",
                protein=protein,
                side=side,
                num_sides=num_sides,
                num_proteins=num_proteins,
            )

        elif edit_type == "UNBIND_NS":
            list_state[index] = 0

        elif edit_type == "BIND_S":
            list_state[index] = self.get_state_id(
                binding_type="BOUND_S",
                protein=protein,
                side=side,
                num_sides=num_sides,
                num_proteins=num_proteins,
            )

        elif edit_type == "SLIDE_RIGHT":
            destination = index + 1

            if destination >= len(list_state):
                raise IndexError(
                    f"Cannot slide right from index {index}."
                )

            list_state[destination] = list_state[index]
            list_state[index] = 0

        elif edit_type == "SLIDE_LEFT":
            destination = index - 1

            if destination < 0:
                raise IndexError(
                    f"Cannot slide left from index {index}."
                )

            list_state[destination] = list_state[index]
            list_state[index] = 0

        else:
            raise ValueError(
                f"Unknown edit type: {edit_type}"
            )

    def close(self) -> None:
        if not self.file_object.closed:
            self.file_object.close()

    def __enter__(self) -> Multi_Trajectory:
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()