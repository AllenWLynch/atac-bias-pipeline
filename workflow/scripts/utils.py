import numpy as np

class lines_grouper:

    def __init__(self, file_object, chunk_size):
        self.file_object = file_object
        self.chunk_size = chunk_size

    def __iter__(self):
        morelines = True
        while morelines:
            group = []
            for i in range(self.chunk_size):
                nextline = self.file_object.readline()
                if nextline == '':
                    morelines = False
                    break
                else:
                    group.append(nextline)
            if len(group) > 0:
                yield group

class Preprocessor:
    
    def __init__(self):
        pass
    
    def fit(self):
        pass
    
    def transform(self, X):
        
        fraglen = X[:,0]
        gc_content = X[:,1]
        biases = X[:,2:]
        
        fraglen_features = np.log2(fraglen)
        bias_features = np.log(biases/(1 - biases))
        
        output = np.hstack([np.vstack([fraglen_features, gc_content]).T, bias_features])

        return output