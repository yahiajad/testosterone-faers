"""Testosterone × sleep apnoea / polycythaemia — FAERS disproportionality, men, stratified by age.
Input: $FAERS_DIR/pq (02_build.py). Output: $FAERS_DIR/results/disproportionality.csv, tto_*.csv; summary printed.
Steps: dedupe (latest version per caseid, drop FDA-deleted cases, collapse exact duplicates on
event_dt/age/sex/country/drug set/PT set) → men only → case/non-case by PT → ROR, PRR+χ², IC (BCPNN).
Sensitivity: PS+SS role; drop lawyer reports (occp_cod LW, 2014–15 testosterone litigation);
drop 2014–2016 reports; active comparator = other drugs used mostly in men (PDE5 inhibitors)."""
import duckdb, os, glob, math, sys
import pandas as pd
from scipy.stats import chi2_contingency
B = os.environ.get("FAERS_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))  # data root (set FAERS_DIR to override)
HERE = os.path.dirname(os.path.abspath(__file__))
P = f"{B}/pq"; OUT = f"{B}/results"; os.makedirs(OUT, exist_ok=True); os.makedirs(f"{B}/tmp", exist_ok=True)
c = duckdb.connect(f"{B}/faers.duckdb")
c.execute(f"set memory_limit='5GB'; set threads=4; set preserve_insertion_order=false; set temp_directory='{B}/tmp'")

TESTO = r"(TESTOST|ANDROGEL|TESTIM|AXIRON|FORTESTA|AVEED|XYOSTED|JATENZO|ANDRODERM|NATESTO|TESTOPEL|VOGELXO|STRIANT|TLANDO|KYZATREX|NEBIDO|SUSTANON)"
TESTO_EXCL = r"(METHYLTESTOST|ESTROGEN|ESTRADIOL|ESTRIOL|PROGESTERONE)"
PDE5 = r"(SILDENAFIL|VIAGRA|TADALAFIL|CIALIS|VARDENAFIL|LEVITRA|AVANAFIL|STENDRA)"
EVENTS = {
    "sleep_apnoea": ["sleep apnoea syndrome", "obstructive sleep apnoea syndrome", "central sleep apnoea syndrome"],
    "polycythaemia": ["polycythaemia", "secondary polycythaemia", "haematocrit increased", "haemoglobin increased",
                      "red blood cell count increased"],
}

def build():
    deleted = set()
    for f in glob.glob(f"{P}/deleted/*.txt"):
        deleted |= {l.strip() for l in open(f) if l.strip()}
    c.execute("create or replace table del as select unnest(?::varchar[]) caseid", [sorted(deleted)])
    # latest version per case
    c.execute(f"""create or replace table demo as
      select * exclude (rn) from (
        select *, row_number() over (partition by caseid order by try_cast(fda_dt as int) desc nulls last,
               try_cast(primaryid as bigint) desc, yq desc, age nulls last, sex nulls last, event_dt nulls last) rn
        from read_parquet('{P}/demo/*.parquet'))
      where rn = 1 and caseid not in (select caseid from del)""")
    c.execute(f"""create or replace table drug as select d.* from read_parquet('{P}/drug/*.parquet') d
                  semi join demo using (primaryid)""")
    c.execute(f"""create or replace table reac as select distinct primaryid, lower(trim(pt)) pt
                  from read_parquet('{P}/reac/*.parquet') r semi join demo using (primaryid)""")
    c.execute(f"""create or replace table ther as select t.* from read_parquet('{P}/ther/*.parquet') t
                  semi join demo using (primaryid)""")
    # collapse exact duplicates across caseids
    # order-independent set fingerprints (sum of distinct hashes) instead of sorted strings: same result, fits in RAM
    c.execute("""create or replace table sig as
        with dg as (select primaryid, sum(distinct hash(upper(coalesce(prod_ai, drugname))))::hugeint ds from drug group by 1),
             rg as (select primaryid, sum(distinct hash(pt))::hugeint ps from reac group by 1)
        select d.primaryid, hash(d.event_dt, d.age, d.sex, d.occr_country, dg.ds, rg.ps) s
        from demo d left join dg using (primaryid) left join rg using (primaryid)""")
    c.execute("""create or replace table base as select d.* from demo d join (
        select min(primaryid) primaryid from sig group by s) k using (primaryid)""")
    c.execute("""create or replace table base as select *,
        case when age_cod='YR' then try_cast(age as double) when age_cod='DEC' then try_cast(age as double)*10
             when age_cod='MON' then try_cast(age as double)/12 when age_cod='WK' then try_cast(age as double)/52
             when age_cod='DY' then try_cast(age as double)/365 end age_y,
        try_cast(substr(coalesce(nullif(init_fda_dt,''), fda_dt),1,4) as int) yr from base""")

def flag(name, pattern, roles, excl=None):
    ex = f"and not regexp_matches(upper(coalesce(prod_ai,'')||' '||drugname), '{excl}')" if excl else ""
    c.execute(f"""create or replace table {name} as select distinct primaryid from drug
        where role_cod in ({",".join(f"'{r}'" for r in roles)})
        and regexp_matches(upper(coalesce(prod_ai,'')||' '||drugname), '{pattern}') {ex}""")

def disprop(a, b, cc, d):
    ror = (a * d) / (b * cc) if b * cc else float("nan")
    se = math.sqrt(1/a + 1/b + 1/cc + 1/d) if min(a, b, cc, d) > 0 else float("nan")
    prr = (a / (a + b)) / (cc / (cc + d)) if cc else float("nan")
    chi = chi2_contingency([[a, b], [cc, d]], correction=True)[0] if min(a + b, cc + d, a + cc, b + d) > 0 else float("nan")
    n = a + b + cc + d; e = (a + b) * (a + cc) / n
    ic = math.log2((a + 0.5) / (e + 0.5))  # shrunk IC (Norén et al.)
    return dict(a=a, b=b, c=cc, d=d, ROR=ror, ROR_lo=math.exp(math.log(ror) - 1.96*se) if a else float("nan"),
                ROR_hi=math.exp(math.log(ror) + 1.96*se) if a else float("nan"), PRR=prr, chi2=chi,
                IC=ic, IC025=ic - 3.3*(a + 0.5)**-0.5 - 2*(a + 0.5)**-1.5,
                signal=bool(a >= 3 and ror and se == se and math.exp(math.log(ror) - 1.96*se) > 1 and ic - 3.3*(a + 0.5)**-0.5 - 2*(a + 0.5)**-1.5 > 0))

def table(exposed, where, label, src="base"):
    rows = []
    for ev, pts in EVENTS.items():
        q = f"""with pop as (select primaryid from {src} where {where}),
          e as (select primaryid from pop semi join {exposed} using (primaryid)),
          o as (select primaryid from pop semi join reac using (primaryid)
                where primaryid in (select primaryid from reac where pt in ({",".join(f"'{p}'" for p in pts)})))
          select (select count(*) from e where primaryid in (select * from o)),
                 (select count(*) from e) - (select count(*) from e where primaryid in (select * from o)),
                 (select count(*) from o) - (select count(*) from e where primaryid in (select * from o)),
                 (select count(*) from pop)"""
        a, b, cc, n = c.execute(q).fetchone(); d = n - a - b - cc
        rows.append(dict(analysis=label, event=ev, **disprop(a, b, cc, d)))
    return rows

if __name__ == "__main__":
    if "--rebuild" in sys.argv or not c.execute("select count(*) from information_schema.tables where table_name='base'").fetchone()[0]:
        build()
    flag("t_ps", TESTO, ["PS"], TESTO_EXCL); flag("t_ss", TESTO, ["PS", "SS"], TESTO_EXCL); flag("pde5", PDE5, ["PS"])
    c.execute(f"""delete from pde5 where primaryid in (
        select primaryid from drug where regexp_matches(upper(drugname), 'REVATIO|ADCIRCA|ALYQ|TADLIQ|LIQREV')
        union select primaryid from read_parquet('{P}/indi/*.parquet') where lower(indi_pt) like '%pulmonary%hypertension%')""")
    flag("t_any", TESTO, ["PS", "SS", "C", "I"], TESTO_EXCL); flag("pde5_any", PDE5, ["PS", "SS", "C", "I"])
    men = "sex='M'"
    R = []
    R += table("t_ps", men, "main: men, PS")
    for lo, hi in [(18, 39), (40, 64), (65, 119)]:
        R += table("t_ps", f"{men} and age_y between {lo} and {hi}.99", f"men {lo}-{hi}, PS")
    R += table("t_ss", men, "sens: PS+SS")
    R += table("t_ps", f"{men} and coalesce(occp_cod,'')<>'LW'", "sens: no lawyer reports")
    R += table("t_ps", f"{men} and yr not between 2014 and 2016", "sens: excl 2014-16")
    # active comparator: testosterone vs PDE5 inhibitors (both male-dominant)
    c.execute("create or replace table comp as select primaryid from t_ps union select primaryid from pde5")
    R += table("t_ps", f"{men} and primaryid in (select primaryid from comp) and primaryid not in (select primaryid from t_any intersect select primaryid from pde5_any)", "active comparator vs PDE5i")
    # sensitivity: skip the exact-duplicate collapse (latest version per case only)
    c.execute("""create or replace table base_nocollapse as select *,
        case when age_cod='YR' then try_cast(age as double) when age_cod='DEC' then try_cast(age as double)*10
             when age_cod='MON' then try_cast(age as double)/12 when age_cod='WK' then try_cast(age as double)/52
             when age_cod='DY' then try_cast(age as double)/365 end age_y,
        try_cast(substr(coalesce(nullif(init_fda_dt,''), fda_dt),1,4) as int) yr from demo""")
    R += table("t_ps", men, "sens: no duplicate collapse", src="base_nocollapse")
    df = pd.DataFrame(R); df.to_csv(f"{OUT}/disproportionality.csv", index=False)
    pd.set_option("display.width", 220); pd.set_option("display.max_columns", 20)
    print(c.execute("select count(*) reports, count(*) filter (where sex='M') men, min(yr), max(yr) from base").df())
    print(df.round(2).to_string(index=False))
    # time to onset (days from therapy start to event), testosterone PS cases with each event
    for ev, pts in EVENTS.items():
        tto = c.execute(f"""select datediff('day', try_strptime(t.start_dt,'%Y%m%d'), try_strptime(b.event_dt,'%Y%m%d')) as tto_days, b.primaryid
            from base b join drug g using (primaryid) join ther t on t.primaryid=b.primaryid and t.dsg_drug_seq=g.drug_seq
            where b.sex='M' and g.role_cod='PS' and b.primaryid in (select * from t_ps)
              and b.primaryid in (select primaryid from reac where pt in ({",".join(f"'{p}'" for p in pts)}))
              and length(t.start_dt)=8 and length(b.event_dt)=8""").df().dropna()
        tto = tto[(tto.tto_days >= 0) & (tto.tto_days <= 7305)]
        tto = tto.groupby("primaryid", as_index=False).tto_days.min()
        tto.to_csv(f"{OUT}/tto_{ev}.csv", index=False)
        if len(tto): print(ev, "time-to-onset n=", len(tto), "median days=", tto.tto_days.median(), "IQR", tto.tto_days.quantile([.25, .75]).tolist())
