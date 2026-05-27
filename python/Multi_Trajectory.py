from __future__ import annotations

import csv
from pathlib import Path


class Multi_Trajectory:
    def __init__(self, csv_name: str) -> None:
        data_dir = Path(__file__).resolve().parent.parent / "data"
        self.file_path = data_dir / csv_name

        if not self.file_path.is_file():
            raise FileNotFoundError(f"CSV file not found: {self.file_path}")

        self.file_object = self.file_path.open(newline="", encoding="utf-8")
        rows = list(csv.reader(self.file_object))
        self.header = rows[0] if rows else []
        self._metadata = self._parse_metadata(self.header)
        self.lines = rows[1:]
        self._current_index = 0

    @staticmethod
    def _parse_metadata(header: list[str]) -> dict[str, str]:
        metadata: dict[str, str] = {}
        for value in header[3:]:
            if "=" not in value:
                continue
            key, parsed_value = value.split("=", 1)
            metadata[key] = parsed_value
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

    def get_runtime(self) -> float:
        return float(self._metadata["runTime"])

    def close(self) -> None:
        if not self.file_object.closed:
            self.file_object.close()

    def applyEdit(self, line: list[str], list_state: list[int]) -> None:
        time, index, edit_type, side = line
        if(edit_type == "BIND_NS" or edit_type == "SWITCH_SIDE" or edit_type == "UNBIND_S"):
            list_state[int(index)] = 2*(int(side)-1) + 1
        elif(edit_type == "UNBIND_NS"):
            list_state[int(index)] = 0
        elif(edit_type == "BIND_S"):
            list_state[int(index)] = 2*(int(side)-1) + 2
        elif(edit_type == "SLIDE_RIGHT"):
            list_state[int(index)+1] = list_state[int(index)]
            list_state[int(index)] = 0
        elif(edit_type == "SLIDE_LEFT"):
            list_state[int(index)-1] = list_state[int(index)]
            list_state[int(index)] = 0

