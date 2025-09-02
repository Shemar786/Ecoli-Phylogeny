#!/usr/bin/env python3
import os, sys, re, shutil, subprocess, glob, csv, tempfile, textwrap

# --------- USER PATHS (edit if needed) ---------
SRC = os.path.expanduser("~/Downloads/E.coli project")  # your genomes folder
WORK = os.path.expanduser("~/Downloads/Ecoli_clean")           # cleaned copies
OUT  = os.path.expanduser("~/Downloads/parsnp_out")            # results
THREADS = "8"
DOCKER_IMAGE = "staphb/parsnp:1.5.6"
# ------------------------------------------------

def fail(msg, code=1):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)

def check_paths():
    if not os.path.isdir(SRC):
        fail(f"Genomes folder not found: {SRC}")
    os.makedirs(WORK, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)

def list_fastas(folder):
    pats = ("*.fasta","*.fa","*.fna","*.fas")
    files = []
    for p in pats:
        files.extend(glob.glob(os.path.join(folder, p)))
    return sorted(files)

def ascii_clean_bytes(in_path, out_path):
    # strip non-ASCII, keep \t \n \r and 0x20-0x7E
    with open(in_path, "rb") as f:
        data = f.read()
    cleaned = bytearray()
    for b in data:
        if b in (9,10,13) or (32 <= b <= 126):
            cleaned.append(b)
    with open(out_path, "wb") as f:
        f.write(cleaned)

def sanitize_fasta(in_path, out_path):
    """
    - ensure ASCII only (already done before this call)
    - normalize headers to safe chars
    - uppercase sequences, keep only A/C/G/T/N
    """
    with open(in_path, "r", errors="ignore") as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                # keep letters, digits, _, . : | - and spaces; collapse spaces to underscore
                h = line.rstrip("\n")
                h = re.sub(r"[^A-Za-z0-9_.:\|\- ]", "_", h)
                h = re.sub(r"\s+", "_", h)
                # avoid empty header
                if h == ">" or h.strip() == "":
                    h = ">unknown"
                fout.write(h + "\n")
            else:
                seq = re.sub(r"[^ACGTNacgtn]", "", line)
                fout.write(seq.upper())

def clean_all_fastas():
    sources = list_fastas(SRC)
    if not sources:
        fail(f"No FASTA files found in {SRC}")
    cleaned = []
    for src in sources:
        base = os.path.basename(src)
        # fix double-dot names like 1..fasta -> 1.fasta
        base = re.sub(r"\.\.f(ast|na|a|as)$", r".fasta", base, flags=re.IGNORECASE)
        # ensure extension is .fasta
        if not re.search(r"\.fasta$", base, re.IGNORECASE):
            base = re.sub(r"\.(fa|fna|fas)$", ".fasta", base, flags=re.IGNORECASE)
        # make sure it's not empty like ".fasta"
        if base.lower() == ".fasta":
            base = "unnamed.fasta"
        tmp_ascii = os.path.join(WORK, f".ascii_{base}")
        dst = os.path.join(WORK, base)
        ascii_clean_bytes(src, tmp_ascii)
        sanitize_fasta(tmp_ascii, dst)
        os.remove(tmp_ascii)
        cleaned.append(dst)
    # pick reference: first cleaned fasta
    cleaned = sorted(cleaned)
    # ensure there is at least one non-empty cleaned file
    kept = []
    for fp in cleaned:
        if os.path.getsize(fp) > 0:
            kept.append(fp)
    if not kept:
        fail("All cleaned FASTAs are empty after sanitization.")
    return kept

def have_docker():
    try:
        subprocess.run(["docker","--version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return True
    except Exception:
        return False

def run_parsnp_docker(ref_path, data_dir, out_dir, threads):
    # mount clean folder and output; run with Linux/amd64 for broad compatibility
    cmd = [
        "docker","run","--rm","--platform=linux/amd64",
        "-v", f"{data_dir}:/data",
        "-v", f"{out_dir}:/out",
        DOCKER_IMAGE,
        "parsnp","-r","/data/"+os.path.basename(ref_path),
        "-d","/data","-o","/out","-p", str(threads)
    ]
    print("\n[INFO] Running Parsnp in Docker:\n" + " ".join(cmd) + "\n")
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        fail("Parsnp (Docker) failed. See messages above.")

def run_parsnp_brew(ref_path, data_dir, out_dir, threads):
    # fallback if user prefers Homebrew-installed parsnp
    parsnp_bin = "/opt/homebrew/bin/parsnp"
    if not os.path.exists(parsnp_bin):
        fail("Homebrew parsnp not found at /opt/homebrew/bin/parsnp and Docker not available.")
    cmd = [
        parsnp_bin, "-r", ref_path, "-d", data_dir, "-o", out_dir, "-p", str(threads)
    ]
    print("\n[INFO] Running Parsnp via Homebrew:\n" + " ".join(f'"{c}"' if " " in c else c for c in cmd) + "\n")
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        fail("Parsnp (Homebrew) failed. See messages above.")

def newick_tip_renamer(newick_in, mapping, newick_out):
    """
    Replace leaf labels exactly matching keys in mapping with values.
    Simple token-wise replace before ':' or delimiters.
    """
    with open(newick_in, "r") as f:
        s = f.read()
    out = []
    tok = ""
    i = 0
    specials = set("(),:; \t\n\r")
    while i < len(s):
        c = s[i]
        if c in "(),;":
            if tok:
                t = tok.strip("'\"")
                out.append(mapping.get(t, tok))
                tok = ""
            out.append(c)
            i += 1
        elif c == ':':
            if tok:
                t = tok.strip("'\"")
                out.append(mapping.get(t, tok))
                tok = ""
            out.append(c)
            i += 1
        elif c.isspace():
            if tok:
                t = tok.strip("'\"")
                out.append(mapping.get(t, tok))
                tok = ""
            out.append(c)
            i += 1
        else:
            tok += c
            i += 1
    with open(newick_out, "w") as f:
        f.write("".join(out))

def build_name_map(clean_files):
    """
    Map each filename stem -> cleaned header name for prettier tree tips.
    """
    mp = {}
    for fp in clean_files:
        stem = os.path.splitext(os.path.basename(fp))[0]
        # read first header
        hdr = "?"
        with open(fp, "r") as f:
            for line in f:
                if line.startswith(">"):
                    hdr = line[1:].strip()
                    break
        # keep it short-ish
        nice = hdr
        # remove trailing words like 'complete_genome' variants for brevity
        nice = re.sub(r"(complete_?genome|chromosome|scaffold|contig|assembly)$","", nice, flags=re.IGNORECASE)
        nice = re.sub(r"__+","_", nice)
        nice = nice.strip("_")
        # fall back to filename if empty
        if not nice:
            nice = stem
        mp[stem] = nice
    return mp

def main():
    print("[INFO] Preparing...")
    check_paths()
    clean_files = clean_all_fastas()
    print(f"[INFO] Cleaned {len(clean_files)} FASTA files into: {WORK}")

    # choose reference (first cleaned file)
    ref = clean_files[0]
    print(f"[INFO] Reference: {ref}")

    # run Parsnp
    if have_docker():
        run_parsnp_docker(ref, WORK, OUT, THREADS)
    else:
        print("[WARN] Docker not found; trying Homebrew parsnp...")
        run_parsnp_brew(ref, WORK, OUT, THREADS)

    # outputs
    tree = os.path.join(OUT, "parsnp.tree")
    xmfa = os.path.join(OUT, "parsnp.xmfa")
    vcf  = os.path.join(OUT, "parsnp.vcf")
    ggr  = os.path.join(OUT, "parsnp.ggr")

    if not os.path.exists(tree):
        print("[WARN] parsnp.tree not found; Parsnp may have filtered too many genomes or failed late.")
    else:
        print(f"[INFO] Tree: {tree}")

    # optional: rename tips to nicer sample names
    if os.path.exists(tree):
        name_map = build_name_map(clean_files)
        # Parsnp typically uses file basenames as tips; build the same keys
        mapping = {}
        for fp in clean_files:
            stem = os.path.splitext(os.path.basename(fp))[0]
            mapping[stem] = name_map.get(stem, stem)

        renamed = os.path.join(OUT, "parsnp_renamed.tree")
        newick_tip_renamer(tree, mapping, renamed)
        print(f"[INFO] Renamed tree (tips prettified): {renamed}")

    print("\n[DONE] Outputs:")
    print(f"  Alignment (XMFA): {xmfa if os.path.exists(xmfa) else '(missing)'}")
    print(f"  SNPs (VCF):       {vcf  if os.path.exists(vcf)  else '(missing)'}")
    print(f"  Tree (Newick):    {tree if os.path.exists(tree) else '(missing)'}")
    print(f"  Tree (renamed):   {os.path.join(OUT,'parsnp_renamed.tree') if os.path.exists(tree) else '(n/a)'}")
    print("\nTip: upload the .tree to iTOL for easy viewing.")

if __name__ == "__main__":
    main()

