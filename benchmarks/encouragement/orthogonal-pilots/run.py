"""Prediction-only budget pilot, disjoint from all evaluation seeds."""
from dataclasses import asdict
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import time
import numpy as np

sys.path.insert(0, '/home/ignacio/longbet-iv-review/benchmarks/encouragement')
from orthogonal_comparison import generate_panel, baseline_features
from longbet._orthogonal_iv import _split_by_assignment
from longbet._iv_nuisance import LongBetIVNuisance, LongBetIVNuisanceConfig, nuisance_validation_score

scenario, seed = sys.argv[1], int(sys.argv[2])
data = generate_panel(scenario, seed)
x = baseline_features(data)
z = data['assignment'].astype(int)
response = np.stack((data['y'], data['d']), axis=-1)[:, data['start']:]
times = data['t'][data['start']:]
fold = _split_by_assignment(z, seed+4189)
train, test = fold == 0, fold == 1
scale = np.maximum(response[train].std(axis=(0, 1)), [1e-8, .05])
record = dict(scenario=scenario, seed=seed, purpose='prediction budget stability only; no IV truth or comparator performance inspected',
    versions={name: importlib.metadata.version(name) for name in ['numpy', 'scipy']}, python=platform.python_version(),
    source_sha256=hashlib.sha256(Path('/home/ignacio/longbet-iv-review/src/longbet/_iv_nuisance.py').read_bytes()).hexdigest(),
    training_accounts=int(train.sum()), test_accounts=int(test.sum()), horizons=len(times), features=x.shape[1], runs=[])
predictions=[]
for burnin, draws in [(100, 200), (200, 400)]:
    config=LongBetIVNuisanceConfig(seed=seed+10000000, burnin=burnin, draws=draws)
    start=time.monotonic()
    model=LongBetIVNuisance(config).fit(x[train], z[train], response[train], times)
    prediction=model.predict(x[test])
    predictions.append(prediction)
    record['runs'].append(dict(config=asdict(config), seconds=time.monotonic()-start,
        metadata=model.metadata_, heldout_observed_arm_normalized_mse=nuisance_validation_score(prediction,z[test],response[test],response[train])))
record['heldout_prediction_rms_change_per_response_training_sd']=np.sqrt(np.mean(((predictions[0]-predictions[1])/scale)**2,axis=(0,1,2))).tolist()
record['heldout_prediction_max_change_per_response_training_sd']=np.max(np.abs((predictions[0]-predictions[1])/scale),axis=(0,1,2)).tolist()
path=Path('/home/ignacio/.cache/longbet-iv-nuisance-pilot')/f'{scenario}_{seed}.json'
path.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'scenario':scenario,'budgets':[{'seconds':r['seconds'],'length_scale':r['metadata']['selected_length_scale'], 'mse':r['heldout_observed_arm_normalized_mse']} for r in record['runs']], 'rms_change':record['heldout_prediction_rms_change_per_response_training_sd']}),flush=True)
