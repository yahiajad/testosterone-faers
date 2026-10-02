"""Convert each FAERS quarterly zip into per-table parquet files ($FAERS_DIR/pq/<table>/<yq>.parquet).
Normalises column names across format changes (BOM, gndr_cod→sex, outc_code→outc_cod), keeps all values as text.
Idempotent: skips quarters already converted."""
import duckdb, zipfile, os, re, glob, tempfile, sys
B = os.environ.get("FAERS_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))  # data root (set FAERS_DIR to override)
HERE = os.path.dirname(os.path.abspath(__file__))
R = f"{B}/raw"; P = f"{B}/pq"
KEEP = {
    "demo": ["primaryid", "caseid", "caseversion", "i_f_code", "event_dt", "init_fda_dt", "fda_dt", "rept_cod",
             "age", "age_cod", "age_grp", "sex", "wt", "wt_cod", "occp_cod", "reporter_country", "occr_country"],
    "drug": ["primaryid", "caseid", "drug_seq", "role_cod", "drugname", "prod_ai", "route", "dose_amt", "dose_unit", "dose_freq"],
    "reac": ["primaryid", "caseid", "pt"],
    "ther": ["primaryid", "caseid", "dsg_drug_seq", "start_dt", "end_dt"],
    "outc": ["primaryid", "caseid", "outc_cod"],
    "indi": ["primaryid", "caseid", "indi_drug_seq", "indi_pt"],
}
RENAME = {"gndr_cod": "sex", "outc_code": "outc_cod"}
con = duckdb.connect()
for z in sorted(glob.glob(f"{R}/faers_ascii_*.zip")):
    yq = re.search(r"(\d{4}q\d)", z).group(1)
    if all(os.path.exists(f"{P}/{t}/{yq}.parquet") for t in KEEP):
        continue
    if not zipfile.is_zipfile(z):  # still downloading
        continue
    with zipfile.ZipFile(z) as zf, tempfile.TemporaryDirectory() as td:
        names = zf.namelist()
        for t, cols in KEEP.items():
            m = [n for n in names if re.search(rf"(^|/){t}\d\dq\d(_new)?\.txt$", n, re.I)]
            if not m:
                print("no", t, yq); continue
            src = zf.extract(m[0], td)
            raw = open(src, "rb").read().replace(b"\xef\xbb\xbf", b"")
            open(src, "wb").write(raw)
            rel = con.sql(f"""select * from read_csv('{src}', delim='$', header=true, all_varchar=true, quote='',
                           escape='', ignore_errors=true, null_padding=true, strict_mode=false)""")
            have = {c.strip().lower(): c for c in rel.columns}
            for a, b in RENAME.items():
                if a in have and b not in have: have[b] = have.pop(a)
            sel = ", ".join(f'"{have[c]}" as {c}' if c in have else f"null::varchar as {c}" for c in cols)
            os.makedirs(f"{P}/{t}", exist_ok=True)
            con.sql(f"copy (select {sel}, '{yq}' as yq from rel) to '{P}/{t}/{yq}.parquet' (format parquet)")
        # FDA-published deletion lists (2019+)
        dl = [n for n in names if re.search(r"delete", n, re.I) and n.lower().endswith(".txt")]
        os.makedirs(f"{P}/deleted", exist_ok=True)
        ids = []
        for n in dl:
            ids += [l.strip() for l in zf.read(n).decode("utf-8", "ignore").splitlines() if l.strip().isdigit()]
        open(f"{P}/deleted/{yq}.txt", "w").write("\n".join(ids))
    print("built", yq, flush=True)
