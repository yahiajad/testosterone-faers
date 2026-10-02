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

# reviewer-requested checks: unknown-age stratum, Mantel-Haenszel age-adjusted ROR, background without any-role testosterone,
# polycythaemia diagnostic vs laboratory-only PTs, narrow polycythaemia definition
a.flag("t_any", a.TESTO, ["PS", "SS", "C", "I"], a.TESTO_EXCL)
def tab(where, pts, extra=""):
    ptl = ",".join(f"'{p}'" for p in pts)
    A, B, C, N = c.execute(f"""with pop as (select primaryid from base where {where} {extra}),
      e as (select primaryid from pop where primaryid in (select * from t_ps)),
      o as (select primaryid from pop where primaryid in (select primaryid from reac where pt in ({ptl})))
      select (select count(*) from e where primaryid in (select * from o)), (select count(*) from e) - (select count(*) from e where primaryid in (select * from o)),
             (select count(*) from o) - (select count(*) from e where primaryid in (select * from o)), (select count(*) from pop)""").fetchone()
    return A, B, C, N - A - B - C
def ror_ci(A, B, C, D):
    r = A * D / (B * C); se = math.sqrt(1/A + 1/B + 1/C + 1/D); return f"{r:.2f} ({math.exp(math.log(r)-1.96*se):.2f}-{math.exp(math.log(r)+1.96*se):.2f}), n={A}"
out.append("\n## Additional checks")
strata = ["age_y < 18", "age_y between 18 and 39.99", "age_y between 40 and 64.99", "age_y between 65 and 119.99", "(age_y is null or age_y >= 120)"]
for ev, pts in a.EVENTS.items():
    na = tab("sex='M' and (age_y is null or age_y >= 120)", pts)
    out.append(f"{ev} age not recorded: ROR {ror_ci(*na)}")
    rows = [tab("sex='M' and " + st, pts) for st in strata]
    T = [(A, B, C, D, A + B + C + D) for A, B, C, D in rows]
    R = sum(A*D/n for A, B, C, D, n in T); S = sum(B*C/n for A, B, C, D, n in T)
    PR = sum((A+D)/n * A*D/n for A, B, C, D, n in T); QS = sum((B+C)/n * B*C/n for A, B, C, D, n in T)
    PSQR = sum((A+D)/n * B*C/n + (B+C)/n * A*D/n for A, B, C, D, n in T)
    mh = R / S; v = PR/(2*R**2) + PSQR/(2*R*S) + QS/(2*S**2)  # Robins-Breslow-Greenland
    out.append(f"{ev} Mantel-Haenszel ROR adjusted for age group (incl. not recorded): {mh:.2f} ({math.exp(math.log(mh)-1.96*math.sqrt(v)):.2f}-{math.exp(math.log(mh)+1.96*math.sqrt(v)):.2f})")
    bg = tab("sex='M'", pts, "and (primaryid in (select * from t_ps) or primaryid not in (select * from t_any))")
    out.append(f"{ev} background excluding testosterone in any role: ROR {ror_ci(*bg)}")
pc = ",".join(f"'{p}'" for p in a.EVENTS["polycythaemia"])
dg, lab = c.execute(f"""with cs as (select primaryid from base where sex='M' and primaryid in (select * from t_ps) and primaryid in (select primaryid from reac where pt in ({pc})))
    select count(*) filter (where primaryid in (select primaryid from reac where pt in ('polycythaemia','secondary polycythaemia'))),
           count(*) filter (where primaryid not in (select primaryid from reac where pt in ('polycythaemia','secondary polycythaemia'))) from cs""").fetchone()
out.append(f"polycythaemia cases with a diagnostic PT: {dg}; laboratory PTs only: {lab}")
nar = tab("sex='M'", ["polycythaemia", "secondary polycythaemia"])
out.append(f"polycythaemia narrow definition (polycythaemia, secondary polycythaemia): ROR {ror_ci(*nar)}")
out.append("sex field, deduplicated reports: " + ", ".join(f"{k}={n}" for k, n in c.execute("select coalesce(sex,'missing'), count(*) from base group by 1 order by 2 desc").fetchall()))
txt = "\n".join(out); open(f"{B}/results/descriptives.txt", "w").write(txt); print(txt)
