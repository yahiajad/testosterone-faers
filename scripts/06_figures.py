"""Figures 1–4 (PNG 300 dpi + PDF) into $FAERS_DIR/figures."""
import os, glob, duckdb, pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
B = os.environ.get("FAERS_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))  # data root (set FAERS_DIR to override)
HERE = os.path.dirname(os.path.abspath(__file__)); F = f"{B}/figures"; os.makedirs(F, exist_ok=True)
plt.rcParams.update({"font.family": "Arial", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
INK, SA, PC = "#222222", "#2b6cb0", "#c05621"
save = lambda fig, n: [fig.savefig(f"{F}/{n}.{x}", dpi=300, bbox_inches="tight") for x in ("png", "pdf")]

# Fig 1: data flow (counts computed from the database built by 03_analysis.py)
import importlib.util
spec = importlib.util.spec_from_file_location("a", f"{HERE}/03_analysis.py"); a = importlib.util.module_from_spec(spec); spec.loader.exec_module(a)
a.flag("t_ps", a.TESTO, ["PS"], a.TESTO_EXCL)
q = lambda sql: a.c.execute(sql).fetchone()[0]
nq = len(glob.glob(f"{a.P}/demo/*.parquet"))
cases = {ev: q(f"""select count(*) from base where sex='M' and primaryid in (select * from t_ps) and primaryid in
              (select primaryid from reac where pt in ({",".join(f"'{p}'" for p in pts)}))""") for ev, pts in a.EVENTS.items()}
n_raw = q(f"select count(*) from read_parquet('{a.P}/demo/*.parquet')")
n_demo = q("select count(*) from demo"); n_base = q("select count(*) from base")
n_men = q("select count(*) from base where sex='M'")
n_t = q("select count(*) from base where sex='M' and primaryid in (select * from t_ps)")
steps = [(f"FAERS quarterly files ({nq} quarters)", f"{n_raw:,} report records"),
         ("Latest version per case; FDA-deleted cases removed", f"{n_demo:,} cases"),
         ("Exact duplicates collapsed (event date, age, sex, country, drug set, PT set)", f"{n_base:,} unique reports"),
         ("Male patients", f"{n_men:,} reports"),
         ("Testosterone as primary suspect", f"{n_t:,} reports"),
         ("Cases: sleep apnoea | polycythaemia", f"{cases['sleep_apnoea']:,} | {cases['polycythaemia']:,}")]
fig, ax = plt.subplots(figsize=(6.2, 6)); ax.axis("off")
for i, (t, n) in enumerate(steps):
    y = 1 - i * 0.18
    ax.add_patch(plt.Rectangle((0.08, y - 0.12), 0.84, 0.11, fill=False, lw=1, ec=INK))
    ax.text(0.5, y - 0.045, t, ha="center", va="center", fontsize=8.5)
    ax.text(0.5, y - 0.09, n, ha="center", va="center", fontsize=9, weight="bold")
    if i < len(steps) - 1: ax.annotate("", (0.5, y - 0.15), (0.5, y - 0.12), arrowprops=dict(arrowstyle="->", color=INK))
ax.set_ylim(-0.1, 1.0); save(fig, "Fig1_flow")

# Fig 2: forest plot
d = pd.read_csv(f"{B}/results/disproportionality.csv")
order = ["main: men, PS", "men 18-39, PS", "men 40-64, PS", "men 65-119, PS", "sens: PS+SS", "sens: no lawyer reports",
         "sens: excl 2014-16", "sens: no duplicate collapse", "active comparator vs PDE5i"]
labels = ["All men (primary analysis)", "  Age 18–39", "  Age 40–64", "  Age ≥65", "Primary or secondary suspect",
          "Lawyer-submitted reports excluded", "First received 2014–2016 excluded", "No exact-duplicate collapse", "Active comparator: PDE5 inhibitors (ED)"]
fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.2), sharey=True)
for ax, ev, col, title in [(axes[0], "sleep_apnoea", SA, "Sleep apnoea"), (axes[1], "polycythaemia", PC, "Polycythaemia")]:
    sub = d[d.event == ev].set_index("analysis").loc[order]
    y = np.arange(len(order))[::-1]
    ax.errorbar(sub.ROR, y, xerr=[sub.ROR - sub.ROR_lo, sub.ROR_hi - sub.ROR], fmt="s", color=col, ms=4, capsize=2, lw=1)
    for yi, (_, r) in zip(y, sub.iterrows()):
        ax.text(1.02, yi, f"{r.ROR:.1f} ({r.ROR_lo:.1f}–{r.ROR_hi:.1f})  n={int(r.a)}", transform=ax.get_yaxis_transform(), va="center", fontsize=7.5)
    ax.axvline(1, color="grey", ls="--", lw=0.8); ax.set_xscale("log"); ax.set_title(title, loc="left", weight="bold")
    ax.set_xlabel("Reporting odds ratio (95% CI), log scale")
axes[0].set_yticks(np.arange(len(order))[::-1]); axes[0].set_yticklabels(labels)
fig.subplots_adjust(wspace=0.9); save(fig, "Fig2_forest")

# Fig 3: cases by FDA receipt year
fig, ax = plt.subplots(figsize=(6.4, 3))
for ev, col, lab in [("sleep_apnoea", SA, "Sleep apnoea"), ("polycythaemia", PC, "Polycythaemia")]:
    pts = ",".join(f"'{p}'" for p in a.EVENTS[ev])
    y = a.c.execute(f"""select yr, count(*) n, count(*) filter (where occp_cod='LW') lw from base where sex='M' and primaryid in (select * from t_ps)
        and primaryid in (select primaryid from reac where pt in ({pts})) and yr>=2012 group by 1 order by 1""").df()
    y = y.set_index("yr").reindex(range(2012, 2027), fill_value=0)
    ax.plot(y.index, y.n, "-o", color=col, ms=3, lw=1.2, label=lab)
    ax.plot(y.index, y.n - y.lw, ":", color=col, lw=1, label=f"{lab}, lawyer reports removed")
ax.axvspan(2013.5, 2016.5, color="grey", alpha=0.12, lw=0, label="2014–2016 (peak lawyer reporting)")
ax.set_xlabel("Year of initial FDA receipt (2026 = Q1–Q2)"); ax.set_ylabel("Cases (men, testosterone PS)"); ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.0, 1.0))
save(fig, "Fig3_by_year")

# Fig 4: time to onset
fig, ax = plt.subplots(figsize=(5, 3))
for ev, col, lab in [("sleep_apnoea", SA, "Sleep apnoea"), ("polycythaemia", PC, "Polycythaemia")]:
    t = pd.read_csv(f"{B}/results/tto_{ev}.csv").tto_days.sort_values()
    ax.step(t / 30.44, np.arange(1, len(t) + 1) / len(t), where="post", color=col, label=f"{lab} (n={len(t)}, median {t.median()/30.44:.1f} mo)")
ax.set_xlabel("Months from testosterone start to event"); ax.set_ylabel("Cumulative proportion of cases"); ax.set_xlim(0, 60)
ax.legend(frameon=False, fontsize=7, loc="lower right"); save(fig, "Fig4_time_to_onset")
print("figures:", sorted(os.listdir(F)))
