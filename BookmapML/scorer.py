import argparse
import json
import os
import traceback
from typing import Dict, List
from tqdm import tqdm
from pathlib import Path

SHOW_PROGRESS_BAR               = True
PRINT_DETAILED_SCORE            = True
T_NUM                           = 3
W                               = [1.5, 1, 0.5]


class FileData:
    def __init__(self):
        self.values = dict()
        self.file_path = ""

    def _validate_probability(self, value):
        return 0 <= value <= 1

    @classmethod
    def read(cls, file_path: str) -> 'FileData':
        file_data = cls()
        file_data.file_path = file_path

        with open(file_path, 'r') as f:
            lines_number = int(f.readline().strip())
            for i in range(lines_number):
                tokens = f.readline().split()

                if len(tokens) != T_NUM + 1:
                    raise Exception(f"Incorrect number of values on line {i + 2} in file {file_path}")
                
                if any([len(x) > 5 for x in tokens[1:]]):
                       raise Exception(f"Some values have more than 3 digits after decimal point on line {i + 2} in file {file_path}")
                
                id = int(tokens[0])
                vals = [float(x) for x in tokens[1:]]
                
                if any([not file_data._validate_probability(x) for x in vals]):
                    raise Exception(f"Invalid probability value on line {i + 2} in file {file_path}")
                
                file_data.values[id] = vals

            if f.read().strip():
                raise Exception(f"File {file_path} contains more lines than required")
            
        return file_data
    
class GroundTruthData(FileData):
    def __init__(self):
        super().__init__()

    def _validate_probability(self, value):
        return value in {-1, 0, 1}
    

class Score:
    def __init__(self, *, score: int = 0, file_path = "<-->"):
        self._score = score
        self._file_path = file_path
        self._scores = [self]

    def accumulate(self, other: 'Score') -> None:
        self._score += other._score
        self._scores += other._scores

    def get_score(self, *, print_detailed_score: bool = True) -> int:
        if print_detailed_score:
            print("scores:")
            print('\n'.join([str(x) for x in self._scores]))
            print('Final score:', self._score)

        return self._score
    
    @classmethod
    def score_for_invalid_file(cls, gt_data: GroundTruthData):
        return Score(file_path = gt_data.file_path)

    @classmethod
    def score(cls, gt_data: GroundTruthData, output_data: FileData) -> 'Score':
        s = 0
        w = 0
        for id in gt_data.values.keys():
            if id in output_data.values:
                output_values = output_data.values[id]
            else:
                output_values = [0.5] * T_NUM                

            gt_values = gt_data.values[id]
            
            for w_t, p_it, g_it in zip(W, output_values, gt_values):
                if g_it == -1:
                    continue
                v = w_t * (2 * g_it - 1) * (2 * p_it - 1)
                s += v
                w += w_t
        
        score = (s + w) / (2 * w)
        score = int(10**7 * max(score, 0))
        return cls(score=score, file_path=output_data.file_path)
    
    def __str__(self):
        return f'{self._file_path}: {self._score}'
    

class Scorer:
    def _get_all_files(self, dir: str) -> Dict[str, str]:
        all_files = {}
        for path, _, files in os.walk(dir):
            for file in files:
                if file in all_files:
                    raise Exception(f"multiple occurences of the same file: {file}")
                all_files[file] = os.path.join(path, file)
        return all_files
    
    def _read_ground_truth_data(self, file_path: str) -> GroundTruthData:
        return GroundTruthData.read(file_path)

    def _read_file_data(self, file_path: str) -> FileData:
        return FileData.read(file_path)

    def _get_output_file_name(self, gt_file_name: str) -> str:
        tokens = gt_file_name.split('.')
        return f"{'.'.join(tokens[:-1])}.out"

    def _read_output_data(self, output_files_paths: Dict[str, str], gt_file_name: str) -> FileData:
        out_file_name = self._get_output_file_name(gt_file_name)
        if not out_file_name in output_files_paths:
            raise Exception(f"file {out_file_name} not found in output directory")
        
        return self._read_file_data(output_files_paths[out_file_name])

    def _get_score_for_absent_file(self, gt_data: GroundTruthData) -> Score:
        return Score.score_for_invalid_file(gt_data)

    def get_score_for_one_file_data(self, gt_data: GroundTruthData, output_data: FileData) -> Score:
        return Score.score(gt_data, output_data)

    def get_score_for_ont_file_paths(self, gt_file_path: str, output_data_path: str) -> Score:
        return self.get_score_for_one_file_data(self._read_file_data(gt_file_path), self._read_file_data(output_data_path))
    
    def _get_score_for_one_file(self, gt_files_paths: Dict[str, str], out_files_paths: Dict[str, str], file_name: str) -> Score:
        gt_data = self._read_ground_truth_data(gt_files_paths[file_name])
        try:
            output_data = self._read_output_data(out_files_paths, file_name)
        except:
            return self._get_score_for_absent_file(gt_data)
        
        return self.get_score_for_one_file_data(gt_data, output_data)

    def get_score(self, gt_dir: str, output_dir: str) -> Score:
        if not os.path.exists(output_dir):
            raise Exception("output directory doesn't exist")
        if not os.path.exists(gt_dir):
            raise Exception("ground truth directory doesn't exist")
        
        gt_files = self._get_all_files(gt_dir)
        out_files = self._get_all_files(output_dir)

        score = Score(file_path="Total")
        
        for file in tqdm(gt_files.keys(), disable=not SHOW_PROGRESS_BAR):
            file_score = self._get_score_for_one_file(gt_files, out_files, file)
            score.accumulate(file_score)

        return score


def main(gt_dir: str, out_dir: str) -> None:
    try:
        score = Scorer().get_score(gt_dir, out_dir)
        print(score.get_score(print_detailed_score=PRINT_DETAILED_SCORE))
    except:
        traceback.print_exc()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--ground-truth-dir', required=True, type=Path, help="The directory with ground truth files (.gt)")
    parser.add_argument('-o', '--output-dir', required=True, type=Path, help="The directory with your output files (.out)")
    args = parser.parse_args()
    main(args.ground_truth_dir, args.output_dir)
