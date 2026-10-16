import copy
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from ercc_audit.core import (InputError, audit_samples, fit_ercc, numeric, linear_fit)
from ercc_audit.cli import main, read_tsv, SAMPLE_FIELDS


def fixture():
    sample=dict(sample_id='SIM-RNA',biological_id='SIM-01',replicate='a',condition='simulation',
        sample_role='sample',rna_ng_ul='12',dna_ng_ul='8',library_input_ul='5',
        recorded_sample_mass_ng='100',ercc_volume_ul='2',ercc_dilution_factor='1000',
        ercc_lot='SIM-LOT',total_effective_bases='10000000')
    blank=dict(sample_id='SIM-NC',biological_id='SIM-NC',replicate='',condition='blank',
        sample_role='negative_control',rna_ng_ul='0',dna_ng_ul='0.166',library_input_ul='14',
        recorded_sample_mass_ng='2.324',ercc_volume_ul='0',ercc_dilution_factor='',
        ercc_lot='',total_effective_bases='')
    standards=[];obs=[]
    for i,(length,copies,covered) in enumerate([(500,2,80),(1000,20,900),(1500,200,1400),(2000,2000,1800)],1):
        rid=f'ERCC-SIM-{i:03d}'
        standards.append(dict(lot_id='SIM-LOT',ercc_id=rid,length_nt=str(length),stock_copies_per_ul=str(copies*500)))
        obs.append(dict(sample_id='SIM-RNA',ercc_id=rid,aligned_bases=str(copies*length/10),
            covered_bases=str(covered),quant_include='1',exclusion_reason=''))
    return [sample,blank],standards,obs


class CalibrationTests(unittest.TestCase):
    def test_actual_mixed_nucleic_acid_mass(self):
        samples,_,_=fixture();row=audit_samples(samples)[0]
        self.assertEqual(row['rna_input_ng'],60)
        self.assertEqual(row['dna_input_ng'],40)
        self.assertEqual(row['calculated_sample_mass_ng'],100)

    def test_nc_cannot_borrow_curve(self):
        args=fixture();details,fits=fit_ercc(*args)
        self.assertEqual(fits[1]['status'],'NO_ERCC')
        self.assertNotIn('y_seq_mean',fits[1])
        self.assertFalse(any(row['sample_id']=='SIM-NC' for row in details))

    def test_known_simulated_curve_and_normalization(self):
        details,fits=fit_ercc(*fixture());fit=fits[0]
        self.assertAlmostEqual(fit['mass_slope'],1)
        self.assertAlmostEqual(fit['mass_r2'],1)
        self.assertAlmostEqual(fit['depth_slope'],1)
        self.assertAlmostEqual(fit['depth_intercept'],1)
        self.assertAlmostEqual(fit['depth_r2'],1)
        self.assertEqual(fit['status'],'EXPLORATORY_FIT_OK')
        self.assertAlmostEqual(details[0]['known_input_copies'],2)
        self.assertAlmostEqual(details[0]['mean_depth'],.2)
        self.assertAlmostEqual(details[0]['breadth'],.16)
        self.assertAlmostEqual(details[0]['y_seq'],details[-1]['y_seq'])

    def test_ercc_mass_includes_excluded_standard(self):
        samples,standards,obs=fixture()
        _,old=fit_ercc(samples,standards,obs)
        obs[-1]['quant_include']='0';obs[-1]['exclusion_reason']='predefined exclusion'
        _,new=fit_ercc(samples,standards,obs)
        self.assertEqual(new[0]['status'],'INSUFFICIENT_LEVELS')
        self.assertEqual(old[0]['all_ercc_mass_ng'],new[0]['all_ercc_mass_ng'])

    def test_missing_observation_is_not_zero(self):
        samples,standards,obs=fixture()
        with self.assertRaisesRegex(InputError,'incomplete observations'):
            fit_ercc(samples,standards,obs[:-1])

    def test_zero_observation_excluded_not_pseudocounted(self):
        samples,standards,obs=fixture()
        obs[0].update(aligned_bases='0',covered_bases='0',quant_include='0',exclusion_reason='not observed')
        details,fits=fit_ercc(samples,standards,obs)
        self.assertEqual(details[0]['y_seq'],0)
        self.assertEqual(fits[0]['status'],'INSUFFICIENT_LEVELS')
        obs[0]['quant_include']='1'
        with self.assertRaisesRegex(InputError,'zero observation'):
            fit_ercc(samples,standards,obs)

    def test_bad_numbers(self):
        for value in ('-1','nan','inf','not a number',''):
            with self.subTest(value=value),self.assertRaises(InputError):numeric(value,'test')
        self.assertIsNone(numeric('','test',optional=True))
        self.assertEqual(numeric('0','test'),0)

    def test_duplicate_keys_and_lot_mismatch(self):
        samples,standards,obs=fixture()
        with self.assertRaisesRegex(InputError,'duplicate sample_id'):
            fit_ercc(samples+[samples[0]],standards,obs)
        with self.assertRaisesRegex(InputError,'duplicate observation'):
            fit_ercc(samples,standards,obs+[obs[0]])
        samples[0]['ercc_lot']='UNKNOWN'
        with self.assertRaisesRegex(InputError,'no matching ERCC standard'):
            fit_ercc(samples,standards,obs)

    def test_inconsistent_coverage_and_denominator(self):
        samples,standards,obs=fixture()
        for field,value in [('covered_bases','500.5'),('covered_bases','501'),('aligned_bases','1')]:
            changed=copy.deepcopy(obs);changed[0][field]=value
            with self.subTest(field=field,value=value),self.assertRaises(InputError):
                fit_ercc(samples,standards,changed)
        samples[0]['total_effective_bases']='100'
        with self.assertRaisesRegex(InputError,'exceed total_effective_bases'):
            fit_ercc(samples,standards,obs)

    def test_missing_total_and_ercc_lot(self):
        samples,standards,obs=fixture();samples[0]['total_effective_bases']=''
        self.assertEqual(audit_samples(samples)[0]['ercc_status'],'MISSING_TOTAL_EFFECTIVE_BASES')
        with self.assertRaisesRegex(InputError,'missing value'):fit_ercc(samples,standards,obs)

    def test_mass_mismatch_exposed(self):
        samples,_,_=fixture();samples[0]['recorded_sample_mass_ng']='60'
        self.assertIn('RECORDED_MASS_MISMATCH',audit_samples(samples)[0]['flags'])

    def test_no_automatic_outlier_deletion(self):
        samples,standards,obs=fixture();obs[1]['aligned_bases']='200000'
        _,fits=fit_ercc(samples,standards,obs)
        self.assertEqual(fits[0]['n_used'],4)
        self.assertEqual(fits[0]['status'],'REVIEW_FIT')

    def test_cli_simulation_and_overwrite_guard(self):
        root=Path(__file__).resolve().parents[1];ex=root/'examples/simulated'
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'run'
            args=['fit-ercc','--samples',str(ex/'samples.tsv'),'--standards',str(ex/'standards.tsv'),
                '--observations',str(ex/'observations.tsv'),'--outdir',str(out)]
            self.assertEqual(main(args),0)
            self.assertTrue((out/'ercc_detail.tsv').is_file())
            manifest=json.loads((out/'run_manifest.json').read_text())
            self.assertEqual(len(manifest['inputs']),3)
            self.assertEqual(main(args),2)

    def test_header_duplicate_and_empty(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'bad.tsv';p.write_text('x\tx\n1\t2\n')
            with self.assertRaisesRegex(InputError,'duplicate column'):read_tsv(p,['x'])
            p.write_text('x\n')
            with self.assertRaisesRegex(InputError,'empty table'):read_tsv(p,['x'])

    def test_constant_fit_refused(self):
        with self.assertRaises(InputError):linear_fit([1,1,1],[1,2,3])


if __name__=='__main__':unittest.main()

