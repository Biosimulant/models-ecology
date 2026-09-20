from pathlib import Path
import hashlib,json,tomllib
import numpy as np
import pytest
from native_fixtures import native_fixture
from src.geci2022_gene_drive import Geci2022GeneDriveModel,GENOTYPES,GAMETES
ROOT=Path(__file__).resolve().parents[3]
FIXTURES=ROOT/'verification'
CASES=json.loads((FIXTURES/'geci-reference-scenarios.json').read_text())


def test_native_fixture_checksums_and_enumeration():
    native=tomllib.loads((FIXTURES/"native-provenance.toml").read_text())
    assert hashlib.sha256((ROOT/"models/core/upstream/MODEL2301120001.jl").read_bytes()).hexdigest()==native["source_sha256"]
    for filename,digest in json.loads((FIXTURES/'fixture-sha256.json').read_text()).items():
        assert hashlib.sha256(native_fixture(filename).read_bytes()).hexdigest()==digest
    assert np.array_equal(np.loadtxt(native_fixture('genotypes.tsv.gz'),dtype=str),np.array(GENOTYPES))
    assert np.array_equal(np.loadtxt(native_fixture('gametes.tsv.gz'),dtype=str),np.array(GAMETES))


@pytest.mark.parametrize('case',CASES,ids=lambda c:c['name'])
def test_every_genotype_and_matrix_against_native_julia(case):
    m=Geci2022GeneDriveModel(**case['kwargs'])
    states=[m._genotype_vector.copy()]
    for generation in range(100):
        m.advance_window(generation,generation+1)
        states.append(m._genotype_vector.copy())
    actual=np.stack(states,axis=1)
    expected=np.loadtxt(native_fixture(case['name']+'-genotypes.tsv.gz'))
    np.testing.assert_allclose(actual,expected,rtol=1e-9,atol=1e-11)
    assert len(m._history)==101
    for key,matrix in m._matrices.items():
        rows=np.loadtxt(native_fixture(f"{case['name']}-{key}_matrix.tsv.gz"),ndmin=2)
        native=np.zeros_like(matrix)
        if matrix.ndim==1:
            assert np.all(rows[:,1]==1)
            native[rows[:,0].astype(int)-1]=rows[:,2]
        else:
            native[rows[:,0].astype(int)-1,rows[:,1].astype(int)-1]=rows[:,2]
        np.testing.assert_allclose(matrix,native,rtol=1e-10,atol=1e-12)


@pytest.mark.parametrize('name,value',[('homing_efficiency',-.1),('editing_efficiency',1.1),('juvenile_survival',1.2),('juvenile_survival',0),('net_reproduction_rate',1),('net_reproduction_rate',.5),('initial_population',0),('release_size',-.1),('recombination_rate',.51),('cas9_cofactor',True),('dominance_editing',float('nan')),('fitness_cost_cas9',float('inf'))])
def test_invalid_constructor_and_input_rejected_atomically(name,value):
    with pytest.raises(ValueError):Geci2022GeneDriveModel(**{name:value})
    m=Geci2022GeneDriveModel()
    with pytest.raises(ValueError):m.set_inputs({'release_size':.2,name:value})
    assert m.release_size==.1


@pytest.mark.parametrize('start,end',[(0,.25),(0,1.5),(0,float('nan')),(0,-1),(1,2)])
def test_invalid_generation_windows_do_not_advance(start,end):
    m=Geci2022GeneDriveModel();before=m._genotype_vector.copy()
    with pytest.raises(ValueError):m.advance_window(start,end)
    assert m._time==0 and not m._history
    np.testing.assert_array_equal(before,m._genotype_vector)


def test_unknown_input_and_parameter_changes_after_start_rejected():
    m=Geci2022GeneDriveModel()
    with pytest.raises(ValueError):m.set_inputs({'missing':1})
    m.advance_window(0,1)
    with pytest.raises(ValueError):m.set_inputs({'homing_efficiency':.8})
    assert m._params['homing_efficiency']==.95
    m.set_inputs({'homing_efficiency':.95})


def test_overrides_apply_before_initial_record_and_reset_repeats():
    m=Geci2022GeneDriveModel();m.set_inputs({'initial_population':2,'release_size':.2})
    m.advance_window(0,5);first=m._genotype_vector.copy()
    assert m._history[0]['total_adults']==pytest.approx(2.4)
    m.reset();m.advance_window(0,5)
    np.testing.assert_array_equal(first,m._genotype_vector)


def test_communication_window_changes_do_not_change_generations():
    fine=Geci2022GeneDriveModel();coarse=Geci2022GeneDriveModel()
    for generation in range(10):fine.advance_window(generation,generation+1)
    coarse.advance_window(0,10)
    assert fine._history==coarse._history


def test_no_release_metric_is_explicitly_change_from_initial_not_treatment_effect():
    m=Geci2022GeneDriveModel(release_size=0);m.advance_window(0,1)
    p=m.get_outputs()['visualisation_payload'].value['payload']
    assert p['history'][1]['total_adults']==pytest.approx(12/7)
    assert p['history'][1]['suppression_ratio']==pytest.approx(-5/7)
    assert 'not a matched no-release treatment effect' in p['interpretation']
    assert m.outputs()['population_state'].emitted_unit=='1'
