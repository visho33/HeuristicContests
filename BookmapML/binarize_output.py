import argparse
import os

#Just binarize the output given a threshold, i.e 
#if p < treshold  -> p = 0
#if p >= treshold ->p = 1
def binarize_file(infile, outfile, threshold):
    
    #Read the file
    with open(infile, 'r') as f:
        lines = f.readlines()

    if not lines:
        return

    try:
        n = int(lines[0].strip())
    except:
        print(f"⚠️ Archivo mal formado: {infile}")
        return

    out_lines = [str(n)]
    for line in lines[1:]:
        parts = line.strip().split()
        if len(parts) != 4:
            continue
        oid, p1, p2, p3 = parts
        #Binarize each probability
        b1 = int(float(p1) >= threshold)
        b2 = int(float(p2) >= threshold)
        b3 = int(float(p3) >= threshold)
        out_lines.append(f"{oid} {b1} {b2} {b3}")

    #Rewrite
    with open(outfile, 'w') as f:
        f.write("\n".join(out_lines) + "\n")

def main():
    
    #Input of the folders
    parser = argparse.ArgumentParser()
    parser.add_argument("--indir", required=True, help="Carpeta con archivos .out (probabilidades)")
    parser.add_argument("--outdir", required=True, help="Carpeta donde guardar los binarizados")
    parser.add_argument("--threshold", type=float, default=0.5, help="Umbral para binarizar (default: 0.5)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    #Binarize each file
    for fname in os.listdir(args.indir):
        if not fname.endswith(".out"):
            continue
        infile = os.path.join(args.indir, fname)
        outfile = os.path.join(args.outdir, fname)
        binarize_file(infile, outfile, args.threshold)

if __name__ == "__main__":
    main()
