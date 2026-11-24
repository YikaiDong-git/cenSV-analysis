#!/usr/bin/env python3
"""frozen_tree_sig — reproducible whole-tree fingerprints for the frozen-result integrity gate."""
import os, sys, json, hashlib, argparse

RECORDED = {
    "step30_phewide":      ("8eba0410e28e98af4074db1ef4309135", "5953812b845b46fdb737319d5872d8c5", 144),
    "step28_phewide":      ("bc16225d3a1591ca150468cf98f54b78", "5900185c61cb40f7edfc08f7f69be3ab", 89),
    "phewide_carryover":   ("cf1e0fdb4b5e35d97bd5e6764f7830a2", "513ddc68f8f38a7c398604d54de09fb0", 10),
    "step30_phewide_mcc1": ("1bf462a3186c83584dbdceaed6f901b6", "b9ae20aca46a354bee699387d0d4da9e", 97),
}
HOME = os.path.expanduser("~")
RV = (f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/"
      f"aou_replication/results_vds")


def _walk(root):
    out = []
    for dirpath, _dirs, files in os.walk(root):
        for fn in files:
            p = os.path.join(dirpath, fn)
            if os.path.islink(p) or not os.path.isfile(p):
                continue
            out.append(os.path.relpath(p, root))
    return sorted(out)


def _md5_file(path, buf=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(buf), b""):
            h.update(chunk)
    return h.hexdigest()


def tree_sig(root):
    rels = _walk(root)
    hc = hashlib.md5(); hm = hashlib.md5()
    for rel in rels:
        p = os.path.join(root, rel)
        st = os.stat(p)
        hc.update(rel.encode()); hc.update(_md5_file(p).encode())
        hm.update(rel.encode()); hm.update(str(st.st_size).encode()); hm.update(str(st.st_mtime_ns).encode())
    return hc.hexdigest(), hm.hexdigest(), len(rels)


def selftest():
    ok = True
    print(f"{'tree':<24}{'n':>5}  {'content_sig':<34}{'meta_sig':<34}verdict")
    for name, (exp_c, exp_m, exp_n) in RECORDED.items():
        root = f"{RV}/{name}"
        if not os.path.isdir(root):
            print(f"{name:<24}{'-':>5}  MISSING TREE"); ok = False; continue
        c, m, n = tree_sig(root)
        vc = "content OK" if c == exp_c else f"CONTENT DRIFT (expected {exp_c})"
        vm = "meta OK" if m == exp_m else "meta CHANGED"
        vn = "" if n == exp_n else f" [n {n} != {exp_n}]"
        if c != exp_c:
            ok = False
        print(f"{name:<24}{n:>5}  {c:<34}{m:<34}{vc}; {vm}{vn}")
    print("\nSELFTEST " + ("PASS — implementation reproduces the recorded fingerprints" if ok
                           else "FAIL — implementation drifted OR a frozen tree changed"))
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("trees", nargs="*", help="tree paths (or bare names under results_vds/)")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--save", default=None, help="write a baseline json")
    ap.add_argument("--verify", default=None, help="compare against a baseline json; exit 1 on any change")
    a = ap.parse_args()

    if a.selftest:
        sys.exit(selftest())

    if a.verify:
        base = json.load(open(a.verify))
        bad = []
        print(f"{'tree':<40}{'n':>6}  verdict")
        for name, rec in sorted(base.items()):
            root = rec["root"]
            if not os.path.isdir(root):
                print(f"{name:<40}{'-':>6}  TREE MISSING"); bad.append(name); continue
            c, m, n = tree_sig(root)
            same_c = c == rec["content_sig"]; same_m = m == rec["meta_sig"]
            if same_c and same_m:
                v = "unchanged (content+meta)"
            elif same_c:
                v = "CONTENT IDENTICAL but meta changed (touched/rewritten, bytes same)"
            else:
                v = "*** CONTENT CHANGED ***"; bad.append(name)
            print(f"{name:<40}{n:>6}  {v}")
        print("\nVERIFY " + ("PASS — no frozen content changed" if not bad
                             else f"FAIL — {len(bad)} tree(s) changed: {bad}"))
        sys.exit(1 if bad else 0)

    out = {}
    print(f"{'tree':<40}{'n':>6}  {'content_sig':<34}meta_sig")
    for t in a.trees:
        root = t if os.path.isdir(t) else f"{RV}/{t}"
        c, m, n = tree_sig(root)
        out[os.path.basename(root.rstrip("/"))] = dict(root=os.path.abspath(root), content_sig=c,
                                                       meta_sig=m, n_files=n)
        print(f"{os.path.basename(root.rstrip('/')):<40}{n:>6}  {c:<34}{m}")
    if a.save:
        json.dump(out, open(a.save, "w"), indent=2, sort_keys=True)
        print(f"\nbaseline -> {a.save}")


if __name__ == "__main__":
    main()
