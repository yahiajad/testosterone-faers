"""Table 1 (characteristics of testosterone-PS cases per event, men) + age-gradient test (ratio of RORs, z-test)
+ reports by year. Reads faers.duckdb and results/disproportionality.csv."""
import duckdb, os, math, pandas as pd
from scipy.stats import norm
import importlib.util
B = os.environ.get("FAERS_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))  # data root (set FAERS_DIR to override)
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("a", f"{HERE}/03_analysis.py"); a = importlib.util.module_from_spec(spec); spec.loader.exec_module(a)
c = a.c
a.flag("t_ps", a.TESTO, ["PS"], a.TESTO_EXCL)
out = []
for ev, pts in a.EVENTS.items():
    ptl = ",".join(f"'{p}'" for p in pts)
    cases = f"select primaryid from base where sex='M' and primaryid in (select * from t_ps) and primaryid in (select primaryid from reac where pt in ({ptl}))"
    q = lambda s: c.execute(s).df()
    n = q(f"select count(*) n from ({cases})").n[0]
    out.append(f"\n## {ev} (n={n})")
    out.append("Age: " + q(f"select count(age_y) known, round(median(age_y),1) median, round(quantile_cont(age_y,.25),1) q1, round(quantile_cont(age_y,.75),1) q3 from base where primaryid in ({cases})").to_string(index=False))
    out.append(q(f"""select case when age_y is null then 'unknown' when age_y<18 then '<18' when age_y<40 then '18-39' when age_y<65 then '40-64' else '65+' end age_grp,
                count(*) n from base where primaryid in ({cases}) group by 1 order by 1""").to_string(index=False))
    out.append(q(f"select coalesce(occp_cod,'unknown') reporter, count(*) n from base where primaryid in ({cases}) group by 1 order by 2 desc").to_string(index=False))
    out.append(q(f"select case when occr_country='US' then 'US' when occr_country is null then 'unknown' else 'non-US' end country, count(*) n from base where primaryid in ({cases}) group by 1 order by 2 desc").to_string(index=False))
    out.append(q(f"""select outc_cod, count(distinct primaryid) n from read_parquet('{a.P}/outc/*.parquet') where primaryid in ({cases}) group by 1 order by 2 desc""").to_string(index=False))
    out.append(q(f"""select upper(coalesce(nullif(route,''),'unknown')) route, count(distinct primaryid) n from drug where role_cod='PS' and primaryid in ({cases})
                 and regexp_matches(upper(coalesce(prod_ai,'')||' '||drugname), '{a.TESTO}') group by 1 order by 2 desc limit 8""").to_string(index=False))
    out.append(q(f"select yr, count(*) n from base where primaryid in ({cases}) group by 1 order by 1").T.to_string(header=False))
out.append("\n## All testosterone-PS reports in men, by age")
out.append(c.execute("""select case when age_y is null then 'unknown' when age_y<18 then '<18' when age_y<40 then '18-39' when age_y<65 then '40-64' else '65+' end age_grp,
    count(*) n from base where sex='M' and primaryid in (select * from t_ps) group by 1 order by 1""").df().to_string(index=False))
# age gradient: ratio of RORs between strata
d = pd.read_csv(f"{B}/results/disproportionality.csv")
def se(r): return math.sqrt(1/r.a + 1/r.b + 1/r.c + 1/r.d)
out.append("\n## Age gradient (ratio of RORs vs 18-39)")
for ev in a.EVENTS:
    ref = d[(d.analysis == "men 18-39, PS") & (d.event == ev)].iloc[0]
    for s in ["men 40-64, PS", "men 65-119, PS"]:
        r = d[(d.analysis == s) & (d.event == ev)].iloc[0]
        lr = math.log(r.ROR / ref.ROR); sd = math.sqrt(se(r)**2 + se(ref)**2)
        out.append(f"{ev:14s} {s:16s} RORR={math.exp(lr):.2f} (95% CI {math.exp(lr-1.96*sd):.2f}-{math.exp(lr+1.96*sd):.2f}) p={2*(1-norm.cdf(abs(lr/sd))):.4f}")
# share of testosterone reports carrying each event, by age (separates numerator from background)
from scipy.stats import fisher_exact
out.append("\n## Event share among testosterone-PS reports and among background, by age")
for ev, pts in a.EVENTS.items():
    ptl = ",".join(f"'{p}'" for p in pts); res = {}
    for lab, lo, hi in [("18-39", 18, 39.99), ("40-64", 40, 64.99), ("65+", 65, 119.99)]:
        r = c.execute(f"""select count(*) filter (where t) nt, count(*) filter (where t and e) et, count(*) filter (where not t) nb, count(*) filter (where e and not t) eb
            from (select primaryid in (select * from t_ps) t, primaryid in (select primaryid from reac where pt in ({ptl})) e
                  from base where sex='M' and age_y between {lo} and {hi})""").fetchone()
        res[lab] = r
        out.append(f"{ev:14s} {lab}: testosterone {r[1]}/{r[0]} = {100*r[1]/r[0]:.2f}% · background {r[3]}/{r[2]} = {100*r[3]/r[2]:.3f}%")
    for lab in ("40-64", "65+"):
        (n0, e0, _, _), (n1, e1, _, _) = res["18-39"], res[lab]
        orr, p = fisher_exact([[e1, n1 - e1], [e0, n0 - e0]])
        out.append(f"   {ev} share {lab} vs 18-39: OR {orr:.2f}, Fisher p={p:.3f}")
txt = "\n".join(out); open(f"{B}/results/descriptives.txt", "w").write(txt); print(txt)
