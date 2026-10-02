"""EBGM / EB05 (DuMouchel 1999 MGPS) for testosterone × target events, men, primary-suspect reports.
Fits a 2-gamma mixture prior on all (PS drug × PT) pairs with N ≥ 1 (zero-truncated likelihood).
Needs faers.duckdb from 03_analysis.py. Signal: EB05 ≥ 2."""
import duckdb, os, numpy as np, pandas as pd
from scipy.optimize import minimize, brentq
from scipy.special import gammaln, digamma
from scipy.stats import gamma as G
import importlib.util
B = os.environ.get("FAERS_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))  # data root (set FAERS_DIR to override)
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("a", f"{HERE}/03_analysis.py"); a = importlib.util.module_from_spec(spec); spec.loader.exec_module(a)
c = duckdb.connect(f"{B}/faers.duckdb")

def counts(where):
    # each target event group is one pseudo-PT, so a report counts once per group (not once per PT)
    grp = " ".join(f"when r.pt in ({','.join(repr(p) for p in pts)}) then '{ev}'" for ev, pts in a.EVENTS.items())
    return c.execute(f"""
      with pop as (select primaryid from base where {where}),
      dr as (select distinct d.primaryid,
               case when regexp_matches(upper(coalesce(prod_ai,'')||' '||drugname), '{a.TESTO}')
                     and not regexp_matches(upper(coalesce(prod_ai,'')||' '||drugname), '{a.TESTO_EXCL}') then 'TESTOSTERONE'
                    else upper(trim(coalesce(nullif(prod_ai,''), drugname))) end drug
             from drug d semi join pop using (primaryid) where role_cod='PS'),
      ev as (select distinct r.primaryid, case {grp} else r.pt end pt from reac r semi join pop using (primaryid)),
      pair as (select drug, pt, count(*) n from dr join ev using (primaryid) group by 1,2),
      nd as (select drug, sum(n) nd from pair group by 1), ne as (select pt, sum(n) ne from pair group by 1),
      tot as (select sum(n) t from pair)
      select drug, pt, n, nd*ne/t e from pair join nd using (drug) join ne using (pt), tot order by drug, pt""").df()

def nb_logpmf(n, alpha, beta, e):
    return gammaln(alpha + n) - gammaln(alpha) - gammaln(n + 1) + alpha*np.log(beta/(beta+e)) + n*np.log(e/(beta+e))

def fit(n, e, k=300_000, seed=1):
    """Prior fitted on a random sample of k pairs (population-level hyperparameters; full-data fit is
    slow and adds nothing), in log space with bounded L-BFGS-B — the naive ratio underflows to 0/0."""
    rng = np.random.default_rng(seed)
    if len(n) > k:
        i = rng.choice(len(n), k, replace=False); n, e = n[i], e[i]
    def comp(al, be):  # zero-truncated NB log-pmf
        l0 = al*np.log(be/(be+e))
        return nb_logpmf(n, al, be, e) - np.log(-np.expm1(l0))
    def nll(x):
        a1, b1, a2, b2 = np.exp(x[:4]); p = 1/(1+np.exp(-x[4]))
        return -np.sum(np.logaddexp(np.log(p) + comp(a1, b1), np.log1p(-p) + comp(a2, b2)))
    best = None
    for x0 in ([np.log(.2), np.log(.1), np.log(2), np.log(4), 0.0], [np.log(1), np.log(1), np.log(5), np.log(5), 1.0], [np.log(.05), np.log(.01), np.log(1.5), np.log(2), -1.0]):
        r = minimize(nll, x0, method="L-BFGS-B", bounds=[(-30, 12)]*4 + [(-12, 12)])
        if best is None or r.fun < best.fun: best = r
    a1, b1, a2, b2 = np.exp(best.x[:4]); p = 1/(1+np.exp(-best.x[4]))
    return a1, b1, a2, b2, p

def post(n, e, th):
    a1, b1, a2, b2, p = th
    l1 = np.log(p) + nb_logpmf(n, a1, b1, e); l2 = np.log(1-p) + nb_logpmf(n, a2, b2, e)
    q = 1/(1+np.exp(l2-l1))
    eblog2 = (q*(digamma(a1+n)-np.log(b1+e)) + (1-q)*(digamma(a2+n)-np.log(b2+e)))/np.log(2)
    cdf = lambda x: q*G.cdf(x, a1+n, scale=1/(b1+e)) + (1-q)*G.cdf(x, a2+n, scale=1/(b2+e))
    eb05 = brentq(lambda x: cdf(x)-0.05, 1e-9, 1e4)
    return 2**eblog2, eb05

if __name__ == "__main__":
    rows = []
    for label, where in [("men, PS", "sex='M'"), ("men 18-39", "sex='M' and age_y between 18 and 39.99"),
                         ("men 40-64", "sex='M' and age_y between 40 and 64.99"), ("men 65+", "sex='M' and age_y between 65 and 119.99")]:
        df = counts(where)
        th = fit(df.n.values.astype(float), df.e.values)
        for ev, pts in a.EVENTS.items():
            sub = df[(df.drug == "TESTOSTERONE") & (df.pt == ev)]
            n, e = float(sub.n.sum()), float(sub.e.sum())  # event group pooled (sum of PT counts)
            if n == 0: rows.append(dict(analysis=label, event=ev, n=0)); continue
            ebgm, eb05 = post(n, e, th)
            rows.append(dict(analysis=label, event=ev, n=n, E=round(e, 2), EBGM=round(ebgm, 2), EB05=round(eb05, 2), signal=eb05 >= 2, at_bound=bool(np.any(np.isclose(np.log(th[:4]), [-30, 12][0]) | np.isclose(np.log(th[:4]), 12))),
                             prior=" ".join(f"{x:.3g}" for x in th), pairs=len(df)))
    out = pd.DataFrame(rows); out.to_csv(f"{B}/results/ebgm.csv", index=False); print(out.to_string(index=False))
