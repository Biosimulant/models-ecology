from pathlib import Path
from types import SimpleNamespace
import importlib.util
import pytest
from biosim import BioWorld,ExecutionPolicy
from src.geci2022_gene_drive import Geci2022GeneDriveModel
p=Path(__file__).resolve().parents[1]/'src/ecology_visualisation.py'
spec=importlib.util.spec_from_file_location('geci_presenter',p)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)


def test_full_graph_presenter_includes_initial_and_final_generation():
    world=BioWorld(communication_step=10)
    core=Geci2022GeneDriveModel()
    v=mod.EcologyVisualisationModel(source_alias='core',mode='gene_drive')
    world.add_biomodule('core',core);world.add_biomodule('visualisation',v)
    world.connect('core.visualisation_payload','visualisation.core_visualisation_payload')
    world.setup();world.run(100)
    assert v.execution_policy is ExecutionPolicy.ONCE_AFTER_RUN
    p=world.get_outputs('core')['visualisation_payload'].value['payload']
    h=p['history'];assert len(h)==101 and h[0]['t']==0 and h[-1]['t']==100
    visuals=v.visualize();plot=visuals[0]['data']
    assert plot['series'][0]['points']==[[p['t'],p['total_adults']] for p in h]
    assert plot['y_unit']=='normalized population'
    assert all(len(row)==len(visuals[2]['data']['columns']) for row in visuals[2]['data']['rows'])
    assert '1 − N(t)/N(0)' in [p['name'] for p in visuals[1]['data']['series']]


@pytest.mark.parametrize('times',[[0,90],[1,100],[0,50,100],[0,float('nan'),100]])
def test_presenter_rejects_incomplete_generation_history(times):
    v=mod.EcologyVisualisationModel(source_alias='core',mode='gene_drive')
    data={'history':[{'t':t,'total_adults':1,'adult_females':.5,'adult_males':.5} for t in times]}
    with pytest.raises(ValueError):v.execute({'core_visualisation_payload':{'payload':data}},context=SimpleNamespace(run_start=0,run_end=100))
