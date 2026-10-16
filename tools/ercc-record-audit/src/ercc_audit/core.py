"""Pure statistical and tabular operations for synthetic ERCC controls."""
import math
import statistics

AVOGADRO = 6.02214076e23
RNA_G_PER_MOL_NT = 340.0  # Approximation from the supplied plan, not lot certification.

class InputError(ValueError):
    pass

def numeric(value, field, *, optional=False, positive=False):
    if value is None or str(value).strip() == "":
        if optional:
            return None
        raise InputError(f"{field}: missing value (blank is not zero)")
    try:
        result = float(value)
    except (ValueError, TypeError):
        raise InputError(f"{field}: expected a number") from None
    if not math.isfinite(result) or result < 0 or (positive and result <= 0):
        raise InputError(f"{field}: expected a finite {'positive' if positive else 'nonnegative'} number")
    return result

def unique(rows, key, label):
    result = {}
    for row in rows:
        k = row.get(key, "").strip()
        if not k or k in result:
            raise InputError(f"{label}: missing or duplicate {key}: {k!r}")
        result[k] = row
    return result

def linear_fit(xs, ys):
    if len(xs) < 3 or len(xs) != len(ys):
        raise InputError("fit needs at least three paired points")
    mx, my = statistics.mean(xs), statistics.mean(ys)
    sxx = sum((x-mx)**2 for x in xs)
    syy = sum((y-my)**2 for y in ys)
    if sxx <= 0 or syy <= 0:
        raise InputError("fit needs variation in both axes")
    slope = sum((x-mx)*(y-my) for x, y in zip(xs, ys)) / sxx
    intercept = my - slope*mx
    residuals = [y-(intercept+slope*x) for x,y in zip(xs,ys)]
    return dict(intercept=intercept, slope=slope,
                r2=max(0.,min(1.,1-sum(r*r for r in residuals)/syy)),
                rmse_log10=math.sqrt(sum(r*r for r in residuals)/len(xs)),
                n=len(xs))

def audit_samples(rows):
    samples = unique(rows, "sample_id", "samples")
    audited=[]
    for sid, row in samples.items():
        for field in ("biological_id", "condition"):
            if not row.get(field, "").strip():
                raise InputError(f"{sid}: missing {field}")
        role = row.get("sample_role", "")
        if role not in ("sample", "negative_control"):
            raise InputError(f"{sid}: sample_role must be sample or negative_control")
        rna=numeric(row.get("rna_ng_ul"), f"{sid}.rna_ng_ul")
        dna=numeric(row.get("dna_ng_ul"), f"{sid}.dna_ng_ul")
        volume=numeric(row.get("library_input_ul"),f"{sid}.library_input_ul",positive=True)
        mass=(rna+dna)*volume
        recorded=numeric(row.get("recorded_sample_mass_ng"),f"{sid}.recorded_sample_mass_ng",optional=True)
        ercc=numeric(row.get("ercc_volume_ul"),f"{sid}.ercc_volume_ul")
        dilution=numeric(row.get("ercc_dilution_factor"),f"{sid}.ercc_dilution_factor",optional=True,positive=True)
        total=numeric(row.get("total_effective_bases"),f"{sid}.total_effective_bases",optional=True,positive=True)
        if mass<=0:
            raise InputError(f"{sid}: zero sample nucleic-acid mass")
        if ercc>0 and (dilution is None or dilution<1):
            raise InputError(f"{sid}: added ERCC requires dilution factor >= 1")
        flags=[]
        if recorded is not None and not math.isclose(mass,recorded,rel_tol=1e-4,abs_tol=1e-6):
            flags.append("RECORDED_MASS_MISMATCH")
        if ercc<=0:
            status="NO_ERCC"
        elif not row.get("ercc_lot", "").strip():
            status="MISSING_ERCC_LOT"
        elif total is None:
            status="MISSING_TOTAL_EFFECTIVE_BASES"
        else:
            status="ERCC_INPUT_READY"
        audited.append(dict(sample_id=sid,biological_id=row['biological_id'],
            replicate=row.get('replicate',''),condition=row['condition'],sample_role=role,
            rna_input_ng=rna*volume,dna_input_ng=dna*volume,
            calculated_sample_mass_ng=mass,recorded_sample_mass_ng=recorded,
            ercc_volume_ul=ercc,ercc_dilution_factor=dilution,
            ercc_status=status,flags=';'.join(flags)))
    return audited

def fit_ercc(sample_rows, standard_rows, observation_rows):
    """Fit ERCC calibration only. Refuse unrecognized or incomplete records."""
    audited=audit_samples(sample_rows)
    samples=unique(sample_rows,'sample_id','samples')
    # Different product lots may legitimately use the same ERCC identifiers.
    standards={}
    for row in standard_rows:
        lot=row.get('lot_id','').strip();rid=row.get('ercc_id','').strip()
        if not lot or not rid.startswith('ERCC-'):
            raise InputError('standards require a lot_id and an ERCC- identifier')
        key=(lot,rid)
        if key in standards:
            raise InputError(f'duplicate ERCC standard: {key}')
        length=numeric(row.get('length_nt'),f'{key}.length_nt',positive=True)
        if length!=int(length):
            raise InputError(f'{key}: length_nt must be an integer')
        concentration=numeric(row.get('stock_copies_per_ul'),f'{key}.stock_copies_per_ul',positive=True)
        standards[key]=dict(length=length,stock=concentration)
    if not standards:
        raise InputError('standards table is empty')
    seen=set();by_sample={sid:[] for sid in samples}
    for row in observation_rows:
        sid=row.get('sample_id','').strip();rid=row.get('ercc_id','').strip()
        if sid not in samples:
            raise InputError(f'observations: unknown sample_id {sid!r}')
        lot=samples[sid].get('ercc_lot','').strip()
        if (lot,rid) not in standards:
            raise InputError(f'{sid}: no matching ERCC standard for {rid!r}, lot {lot!r}')
        key=(sid,rid)
        if key in seen:
            raise InputError(f'duplicate observation {key}')
        seen.add(key)
        bases=numeric(row.get('aligned_bases'),f'{key}.aligned_bases')
        covered=numeric(row.get('covered_bases'),f'{key}.covered_bases')
        if covered!=int(covered) or covered>standards[(lot,rid)]['length']:
            raise InputError(f'{key}: covered_bases must be an integer within reference length')
        if bases<covered or ((bases==0)!=(covered==0)):
            raise InputError(f'{key}: aligned_bases and covered_bases are inconsistent')
        included=row.get('quant_include','').strip()
        if included not in ('0','1'):
            raise InputError(f'{key}: quant_include must be 0 or 1')
        reason=row.get('exclusion_reason','').strip()
        if included=='0' and not reason:
            raise InputError(f'{key}: excluded observations require a reason')
        by_sample[sid].append(dict(ercc_id=rid,aligned_bases=bases,
            covered_bases=covered,include=included=='1',exclusion_reason=reason))
    detail=[];fits=[]
    for a in audited:
        sid=a['sample_id'];row=samples[sid]
        if a['ercc_volume_ul']==0:
            if by_sample[sid]:
                raise InputError(f'{sid}: observations provided for a sample without ERCC input')
            fits.append(dict(sample_id=sid,status='NO_ERCC',n_used=0))
            continue
        lot=row.get('ercc_lot','').strip()
        if not lot:
            raise InputError(f'{sid}: missing ercc_lot')
        lot_std={rid:val for (l,rid),val in standards.items() if l==lot}
        if not lot_std:
            raise InputError(f'{sid}: no standards for lot {lot!r}')
        obs={o['ercc_id']:o for o in by_sample[sid]}
        if set(obs)!=set(lot_std):
            missing=sorted(set(lot_std)-set(obs))
            raise InputError(f'{sid}: incomplete observations; explicit zero rows required; missing {missing}')
        total=numeric(row.get('total_effective_bases'),f'{sid}.total_effective_bases',positive=True)
        sum_bases=sum(o['aligned_bases'] for o in obs.values())
        if sum_bases>total:
            raise InputError(f'{sid}: summed ERCC aligned bases exceed total_effective_bases')
        vlib=numeric(row['library_input_ul'],'library_input_ul',positive=True)
        scale=a['ercc_volume_ul']/a['ercc_dilution_factor']
        # Include ALL supplied standards in total ERCC mass, including excluded ones.
        ercc_mass=sum(s['stock']*scale*s['length']*RNA_G_PER_MOL_NT/AVOGADRO*1e9 for s in lot_std.values())
        total_mass=a['calculated_sample_mass_ng']+ercc_mass
        local=[]
        for rid,s in lot_std.items():
            o=obs[rid];copies=s['stock']*scale
            mass=copies*s['length']*RNA_G_PER_MOL_NT/AVOGADRO*1e9
            d=dict(sample_id=sid,ercc_id=rid,lot_id=lot,length_nt=int(s['length']),
                known_input_copies=copies,known_input_mass_ng=mass,
                known_equivalent_copies_per_ul=copies/vlib,
                aligned_bases=o['aligned_bases'],mean_depth=o['aligned_bases']/s['length'],
                breadth=o['covered_bases']/s['length'],mass_fraction=mass/total_mass,
                base_fraction=o['aligned_bases']/total,
                quant_include=int(o['include']),exclusion_reason=o['exclusion_reason'])
            d['y_seq']=(d['base_fraction']/d['mass_fraction']) if d['aligned_bases']>0 else 0.
            # This is a dimensionally normalized depth, not measured concentration.
            d['normalized_depth_per_ul']=d['mean_depth']/vlib
            if o['include'] and d['aligned_bases']<=0:
                raise InputError(f'{sid}/{rid}: cannot include a zero observation in a log fit')
            local.append(d)
        used=[d for d in local if d['quant_include']]
        levels=len(set(d['known_input_copies'] for d in used))
        result=dict(sample_id=sid,status='INSUFFICIENT_LEVELS',n_used=len(used),
            concentration_levels=levels,sample_mass_ng=a['calculated_sample_mass_ng'],
            all_ercc_mass_ng=ercc_mass,total_input_mass_ng=total_mass)
        if levels>=4:
            mass_fit=linear_fit([math.log10(d['mass_fraction']) for d in used],
                                [math.log10(d['base_fraction']) for d in used])
            depth_fit=linear_fit([math.log10(d['normalized_depth_per_ul']) for d in used],
                                [math.log10(d['known_equivalent_copies_per_ul']) for d in used])
            for prefix, fit in [('mass',mass_fit),('depth',depth_fit)]:
                for key,val in fit.items():result[f'{prefix}_{key}']=val
            result['y_seq_mean']=statistics.mean(d['y_seq'] for d in used)
            result['y_seq_median']=statistics.median(d['y_seq'] for d in used)
            result['depth_observed_min']=min(d['normalized_depth_per_ul'] for d in used)
            result['depth_observed_max']=max(d['normalized_depth_per_ul'] for d in used)
            # Plan-defined exploratory checks; no auto-exclusion and no validated LoQ.
            ok=(mass_fit['r2']>=.95 and .8<=mass_fit['slope']<=1.2 and
                depth_fit['r2']>=.95 and .8<=depth_fit['slope']<=1.2)
            result['status']='EXPLORATORY_FIT_OK' if ok else 'REVIEW_FIT'
        detail.extend(local);fits.append(result)
    return detail,fits

