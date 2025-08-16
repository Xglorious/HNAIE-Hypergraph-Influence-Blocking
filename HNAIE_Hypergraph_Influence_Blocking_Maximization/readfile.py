
import threading

class DataLoader:
    def __init__(self, file_path):
        self.file_path = file_path
        self.output_file = None
        self.lock = threading.Lock() 

    def read_hypergraph_from_txt(self):
        with open(self.file_path, 'r') as file:
            lines = file.readlines()

        vertices = set()
        hyperedges = {}

        for i, line in enumerate(lines):
            parts = line.strip().split()
            hyperedge = set(parts)
            hyperedges[f'e{i}'] = hyperedge
            vertices.update(hyperedge)

        return vertices, hyperedges

    def open_output_file(self, write_file_path):
        with self.lock:
            if self.output_file is not None:
                self.output_file.close()
            self.output_file = open(write_file_path, 'w', encoding='utf-8')

    def print_to_file(self, *args):
        with self.lock:
            if self.output_file is None:
                raise ValueError("Output file is not opened. Please call open_output_file() first.")
            print(*args, file=self.output_file)
            self.output_file.flush()

    def close_output_file(self):
        with self.lock:
            if self.output_file:
                self.output_file.close()
                self.output_file = None