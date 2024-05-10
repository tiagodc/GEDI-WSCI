import os, dill, sys, argparse, re
import numpy as np
import pandas as pd
import xgboost as xgb
from crepes import WrapRegressor
from crepes.extras import binning
import pynndescent as pnn

MODELS_ROOT = './models'

def getCmdArgs():
    p = argparse.ArgumentParser(description = "Apply Waveform Structural Complexity (WSCI) model to input RH metrics")
    
    p.add_argument("-i", "--input", dest="input", required=True, type=str, help="input RH metrics file [cm]")
    p.add_argument("-o", "--output", dest="output", required=True, type=str, help="output file")
    p.add_argument("-p", "--pft", dest="pft", required=False, type=str, default=None, help='PFT code or column name')
    p.add_argument("-r", "--rh-cols", dest="rh_cols", required=False, type=str, default=None, help='string pattern to select RH columns')
    p.add_argument("-x", "--index", dest="index", required=False, type=str, default=[], nargs='+', help='columns to copy from input file into the output')
    
    cmdargs = p.parse_args()
    return cmdargs


class WSCI:
    def __init__(self, regressor, neighbors=25, X_cal=None, y_cal=None):
        self.regressor = regressor
        self.neighbors = neighbors        
        if X_cal is not None and y_cal is not None:
            self.store_calibration_data(X_cal, y_cal)

    def prepare(self, X_train, y_train, jobs=10):
        self.nn = pnn.NNDescent(X_train, n_neighbors=25, n_jobs=10)
        self.nn.prepare()
        self.residuals = y_train - self.regressor.predict(X_train)

    def apply(self, X):
        knn = self.nn.query(X, self.neighbors)
        ids = knn[0]
        sigmas = np.array([np.mean(np.abs(self.residuals[i])) for i in ids])
        return sigmas 
            
    def store_calibration_data(self, X_cal, y_cal):
        self.X_cal = X_cal
        self.y_cal = y_cal

    def fit(self, alpha=.05, n_bins=0, use_sigmas=False):
        self.bin_thresholds = None  
        self.use_sigmas = use_sigmas
        kw = {}
        
        if n_bins > 0:
            bins_cal, bin_thresholds = binning(self.regressor.predict(self.X_cal), bins=n_bins)
            self.bin_thresholds = bin_thresholds
            kw['bins'] = bins_cal
        if use_sigmas:
            kw['sigmas'] = self.apply(X_cal)

        mod = WrapRegressor(self.regressor)
        self.conformal_regressor = mod.calibrate(self.X_cal, self.y_cal, **kw)
        return self.conformal_regressor.evaluate(self.X_cal, self.y_cal, **kw)

    def predict(self, X, confidence=0.95):
        preds = self.regressor.predict(X)
        
        if confidence == 0:
            return preds
        
        kw = {}
        if self.bin_thresholds is not None:
            from crepes.extras import binning
            kw['bins'] = binning(preds, self.bin_thresholds)
        if self.use_sigmas:
            kw['sigmas'] = self.apply(X)

        cf_preds = self.conformal_regressor.predict_int(X, confidence=confidence, y_min=0, **kw)
        return preds, cf_preds, kw

    def clear_data(self):
        del self.X_cal
        del self.y_cal

    def save(self, path):
        with open(path, 'wb') as model_file:
            dill.dump(self, model_file)
            
def load_predictor(path):
    with open(path, 'rb') as f:
        predictor = dill.load(f)
        return predictor
            
def predict_pft(xdf, rh_cols, pft_col='pft_class', idx_cols=[]):
    pft = xdf[pft_col].iloc[0]
    if pft in [5,6,11]: pft = 5611
    
    mod_suffix =  f"pft_{pft}" if pft in [1,2,4,5611] else 'global'    
    mod_path = f"{MODELS_ROOT}/wsci_model_{mod_suffix}.pkl" 
    mod = load_predictor(mod_path)    
    
    wsci, wsci_err, kw = mod.predict(xdf[rh_cols].to_numpy(), confidence=0.95)    
    xdf = xdf.assign(wsci = wsci, wsci_pi = wsci_err[:,1] - wsci)
    wsci_cols = ['wsci', 'wsci_pi']
        
    for i in ['XY','Z']: 
        mod_path = f"{MODELS_ROOT}/wsci_model_{i}_{mod_suffix}.pkl"         
        mod = load_predictor(mod_path)    
        wsci, wsci_err, kw = mod.predict(xdf[rh_cols].to_numpy(), confidence=0.95)    
        icols = {f"wsci_{i.lower()}": wsci, f"wsci_{i.lower()}_pi": wsci_err[:,1] - wsci}
        wsci_cols += icols.keys()
        xdf = xdf.assign(**icols)            
    
    return xdf[idx_cols + wsci_cols]

def apply_rh_model(df, rh_cols, pft_col='pft_class', idx_cols=[]):            
    wsci_pft = df.groupby(pft_col, group_keys=False).apply(predict_pft, rh_cols=rh_cols, idx_cols=idx_cols)
    return wsci_pft

if __name__ == '__main__':    
    args = getCmdArgs()
       
    print("## -- reading rh metrics")
    if args.input.endswith('.parquet') or args.input.endswith('.pq') or args.input.endswith('.parq'):
        rhdf = pd.read_parquet(args.input)
    elif args.input.endswith('.csv'):
        rhdf = pd.read_csv(args.input)
    else:
        rhdf = pd.read_table(args.input, delim_whitespace=True, header=None, index_col=None)
    
    if args.rh_cols is not None:
        rex = re.compile(args.rh_cols)
        rh_cols = [i for i in rhdf.columns if rex.match(i) is not None]
    else:
        rh_cols = rhdf.columns[:101]
        
    n = len(rh_cols)
    if n < 101:
        sys.exit('## -- input RH metrics incomplete, needs 101 RH columns')
    
    print("## -- extracting PFTs")
    pft_col = 'pft_class'    
    if args.pft is None:
        rhdf[pft_col] = -1
    elif args.pft in rhdf.columns:
        pft_col = args.pft
    else:
        rhdf[pft_col] = int(args.pft)

    print("## -- calculating WSCI")
    wsci_df = apply_rh_model(rhdf, rh_cols, pft_col, args.index)
    wsci_df = wsci_df.loc[rhdf.index]

    print("## -- exporting results")
    fmt = args.output.split('.')[-1]
    if fmt in ['parquet','pq''parq']:
        wsci_df.to_parquet(args.output)
    else:
        wsci_df.to_csv(args.output, index=False)
    
    sys.exit('## -- DONE')