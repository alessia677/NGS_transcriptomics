docker run -i -t -v /Users/ferdinandosquitieri/Desktop/transcriptomics/tutorial:/tutorial -w /tutorial ceciliaklein/teaching:uvic

export PATH=$PATH:/tutorial/teaching-utils/;

1. Perform differential expression analysis between brain and liver using the EdgeR. Present results using a heatmap with hierarchical clustering in rows and columns and colored classification of differentially expressed genes (DEGs), i.e. overexpressed in
brain versus liver and the other way around.

cd analysis

edgeR.analysis.R --input_matrix ../quantifications/encode.mouse.gene.expected_count.idr_NA.tsv \
                 --metadata /tutorial/data/gene.quantifications.index.tsv \
                 --fields tissue \
                 --coefficient 3 \
                 --output brain_X_liver
                 

1. DEGs IN BRAIN (118)
awk '$NF<0.01 && $2<-10{print $1"\tover_brain_X_liver"}' edgeR.cpm1.n2.brain_X_liver.tsv > edgeR.over_brain_X_liver.txt

2. DEGs IN LIVER (107)
awk '$NF<0.01 && $2>10 {print $1"\tover_liver_X_brain"}' edgeR.cpm1.n2.brain_X_liver.tsv > edgeR.over_liver_X_brain.txt

- HEATMAP:

awk '$3=="gene"{ match($0, /gene_id "([^"]+).+gene_type "([^"]+)/, var); print var[1],var[2] }' OFS="\t" /tutorial/refs/gencode.vM4.gtf \
| join.py --file1 stdin \
          --file2 <(cat edgeR.over*.txt) \
| sed '1igene\tedgeR\tgene_type' > gene.edgeR.2.tsv


cut -f1 gene.edgeR.2.tsv \
| tail -n+2 \
| selectMatrixRows.sh - ../quantifications/encode.mouse.gene.TPM.idr_NA.tsv \
| ggheatmap.R --width 5 \
              --height 8 \
              --col_metadata /tutorial/data/gene.quantifications.index.tsv \
              --colSide_by tissue \
              --col_labels labExpId \
              --row_metadata gene.edgeR.2.tsv \
              --merge_row_mdata_on gene \
              --rowSide_by edgeR,gene_type \
              --row_labels none \
              --log \
              --pseudocount 0.1 \
              --col_dendro \
              --row_dendro \
              --matrix_palette /tutorial/palettes/palDiverging.txt \
              --colSide_palette /tutorial/palettes/palTissue.txt \
              --output heatmap.brain_X_liver.pdf

2. Perform gene ontology enrichment analysis of the two sets of DEGs using the command line wrapper of GOstats R package for biological processes. Plot results using any
graphical representation and discuss results.

awk '{split($10,a,/\"|\./); print a[2]}' /tutorial/refs/gencode.vM4.gtf | sort -u > universe.txt

1.  BP brain X liver
awk '{split($1,a,"."); print a[1]}' edgeR.over_brain_X_liver.txt \
| GO_enrichment.R --universe universe.txt \
                  --genes stdin \
                  --categ BP \
                  --output edgeR.over_brain_X_liver \
                  --species mouse

awk 'NR==1{$1="% "$1}{print $1,$2}' edgeR.over_brain_X_liver.BP.tsv

2. BP liver X brain

awk '{split($1,a,"."); print a[1]}' edgeR.over_liver_X_brain.txt \
| GO_enrichment.R --universe universe.txt \
                  --genes stdin \
                  --categ BP \
                  --output edgeR.over_liver_X_brain \
                  --species mouse

awk 'NR==1{$1="% "$1}{print $1,$2}' edgeR.over_liver_X_brain.BP.tsv

3. Analyze differential splicing using SUPPA between brain and liver for skipping exon, intron retention, mutually exclusive exon and alternative first exon. Plot top results using
heatmaps. Different thresholds may be chosen for each event type.

cd ../splicing

- LIST OF TRANSCRIPT ID -
awk '$3=="transcript" && $0~/gene_type "protein_coding"/{ match($0, /transcript_id "([^"]+)/, id); print id[1] }' /tutorial/refs/gencode.vM4.gtf |sort -u > protein_coding_transcript_IDs.txt

- Genome annotation restricted to exon features and filtered by transcript type -
cat /tutorial/refs/gencode.vM4.gtf |awk '$3=="exon"' |grep -Ff protein_coding_transcript_IDs.txt > exon-annot.gtf

- Filter transcript TPM matrix -
selectMatrixRows.sh protein_coding_transcript_IDs.txt /tutorial/quantifications/encode.mouse.transcript.TPM.idr_NA.tsv > pc-tx.tsv

- Individual transcript expression matrices -
for tissue in Brain Heart Liver; do
    selectMatrixColumns.sh PRNAembryo${tissue}1:PRNAembryo${tissue}2 pc-tx.tsv > expr.${tissue}.tsv
done

- generating SE, MX, FL and RI events -
suppa.py generateEvents -i exon-annot.gtf -e SE MX RI FL -o localEvents -f ioe

- Compute percent spliced in index (PSI) values for local events -

1. EVENT SE:

event=SE; suppa.py psiPerEvent --total-filter 10 --ioe-file localEvents_${event}_strict.ioe --expression-file pc-tx.tsv -o PSI-${event}

event=SE; for tissue in Brain Heart Liver ;do selectMatrixColumns.sh PRNAembryo${tissue}1:PRNAembryo${tissue}2 PSI-${event}.psi > ${tissue}.${event}.psi;done

2. EVENT RI:
event=RI; suppa.py psiPerEvent --total-filter 10 --ioe-file localEvents_${event}_strict.ioe --expression-file pc-tx.tsv -o PSI-${event}

event=RI; for tissue in Brain Liver ;do selectMatrixColumns.sh PRNAembryo${tissue}1:PRNAembryo${tissue}2 PSI-${event}.psi > ${tissue}.${event}.psi;done

3. EVENT MX:

event=MX; suppa.py psiPerEvent --total-filter 10 --ioe-file localEvents_${event}_strict.ioe --expression-file pc-tx.tsv -o PSI-${event}

event=MX; for tissue in Brain Liver ;do selectMatrixColumns.sh PRNAembryo${tissue}1:PRNAembryo${tissue}2 PSI-${event}.psi > ${tissue}.${event}.psi;done

4. EVENT AF:

event=AF; suppa.py psiPerEvent --total-filter 10 --ioe-file localEvents_${event}_strict.ioe --expression-file pc-tx.tsv -o PSI-${event}

event=AF; for tissue in Brain Liver ;do selectMatrixColumns.sh PRNAembryo${tissue}1:PRNAembryo${tissue}2 PSI-${event}.psi > ${tissue}.${event}.psi;done

- Differential splicing analysis for local events -

event=SE; suppa.py diffSplice --method empirical --input localEvents_${event}_strict.ioe --psi Brain.${event}.psi Heart.${event}.psi Liver.${event}.psi --tpm expr.Brain.tsv expr.Heart.tsv expr.Liver.tsv -c -gc -o DS.${event}

event=SE; awk 'BEGIN{FS=OFS="\t"}NR==1{print }NR>1 && $2!="nan" && ($2>0.5 || $2<-0.5) && $3<0.05{print}' DS.${event}.dpsi

event=RI; suppa.py diffSplice --method empirical --input localEvents_${event}_strict.ioe --psi Brain.${event}.psi Liver.${event}.psi --tpm expr.Brain.tsv expr.Liver.tsv -c -gc -o DS.${event}

event=RI; awk 'BEGIN{FS=OFS="\t"}NR==1{print }NR>1 && $2!="nan" && ($2>0.3 || $2<-0.3) && $3<0.05{print}' DS.${event}.dpsi

event=MX; suppa.py diffSplice --method empirical --input localEvents_${event}_strict.ioe --psi Brain.${event}.psi Liver.${event}.psi --tpm expr.Brain.tsv expr.Liver.tsv -c -gc -o DS.${event}

event=MX; awk 'BEGIN{FS=OFS="\t"}NR==1{print }NR>1 && $2!="nan" && ($2>0.3 || $2<-0.3) && $3<0.05{print}' DS.${event}.dpsi

event=AF; suppa.py diffSplice --method empirical --input localEvents_${event}_strict.ioe --psi Brain.${event}.psi Liver.${event}.psi --tpm expr.Brain.tsv expr.Liver.tsv -c -gc -o DS.${event}

event=AF; awk 'BEGIN{FS=OFS="\t"}NR==1{print }NR>1 && $2!="nan" && ($2>0.5 || $2<-0.5) && $3<0.05{print}' DS.${event}.dpsi

- HEATMAP -

Heatmap:

SE: 

# prepare input for heatmap
event=SE; awk 'BEGIN{FS=OFS="\t"}NR>1 && $2!="nan" && ($2>0.5 || $2<-0.5) && $3<0.05{print}' DS.${event}.dpsi|cut -f1 > top-examples-SE.txt
selectMatrixRows.sh top-examples-SE.txt DS.SE.psivec > matrix.top-examples-SE.tsv

# heatmap SE
ggheatmap.R -i matrix.top-examples-SE.tsv -o heatmap_top-examples-SE.pdf --matrix_palette /tutorial/palettes/palSequential.txt --row_dendro  --matrix_fill_limits "0,1" -B 8

RI: 

# prepare input for heatmap
event=RI; awk 'BEGIN{FS=OFS="\t"}NR>1 && $2!="nan" && ($2>0.3 || $2<-0.3) && $3<0.05{print}' DS.${event}.dpsi|cut -f1 > top-examples-RI.txt
selectMatrixRows.sh top-examples-RI.txt DS.RI.psivec > matrix.top-examples-RI.tsv

# heatmap SE
ggheatmap.R -i matrix.top-examples-RI.tsv -o heatmap_top-examples-RI.pdf --matrix_palette /tutorial/palettes/palSequential.txt --row_dendro  --matrix_fill_limits "0,1" -B 8


MX: 

# prepare input for heatmap
event=MX; awk 'BEGIN{FS=OFS="\t"}NR>1 && $2!="nan" && ($2>0.3 || $2<-0.3) && $3<0.05{print}' DS.${event}.dpsi|cut -f1 > top-examples-MX.txt
selectMatrixRows.sh top-examples-MX.txt DS.MX.psivec > matrix.top-examples-MX.tsv

# heatmap SE
ggheatmap.R -i matrix.top-examples-MX.tsv -o heatmap_top-examples-MX.pdf --matrix_palette /tutorial/palettes/palSequential.txt --row_dendro  --matrix_fill_limits "0,1" -B 8

AF: 

# prepare input for heatmap
event=AF; awk 'BEGIN{FS=OFS="\t"}NR>1 && $2!="nan" && ($2>0.6 || $2<-0.6) && $3<0.05{print}' DS.${event}.dpsi|cut -f1 > top-examples-AF.txt
selectMatrixRows.sh top-examples-AF.txt DS.AF.psivec > matrix.top-examples-AF.tsv

# heatmap SE
ggheatmap.R -i matrix.top-examples-AF.tsv -o heatmap_top-examples-AF.pdf --matrix_palette /tutorial/palettes/palSequential.txt --row_dendro  --matrix_fill_limits "0,1" -B 8

4. Find H3K4me3 peaks shared by brain and liver and the ones exclusively found in each tissue using the narrow peaks found in /tutorial/results using bedtools intersect. Show results using a bar plot colored by the color code used during the hands-on. Palette is
availables at /tutorial/palettes/palTissue.txt. Any color may be chosen for shared peaks.

mkdir ../chip-analysis

cd chip-analysis

bedtools intersect -a /tutorial/results/CHIPembryoBrain.narrowPeak -b /tutorial/results/CHIPembryoLiver.narrowPeak -wao|awk 'BEGIN{FS=OFS="\t"; print "Brain_coordinates","Brain_peak","Liver_coordinates","Heart_peak","intersection"}$NF!=0{printf ("%20s\t%15s\t%15s\t%10s\t%20s\t%15s\t%15s\t%10s\t\n", $1,$2,$3,$4,$11,$12,$13,$14)}' > common-peaks-Brain-liver.tsv

bedtools intersect -a /tutorial/results/CHIPembryoBrain.narrowPeak -b /tutorial/results/CHIPembryoLiver.narrowPeak -wao|awk 'BEGIN{OFS="\t"; print "Brain_coordinates","Brain_peak","Liver_coordinates","Liver_peak","intersection"}$NF!=0{print $1,$2,$3,$4,$11,$12,$13,$14}' > Common-peaks-Brain-liver.bed

bedtools intersect -a /tutorial/results/CHIPembryoBrain.narrowPeak -b /tutorial/results/CHIPembryoLiver.narrowPeak -wao|awk 'BEGIN{FS=OFS="\t"; print "Brain_coordinates","Brain_peak"}$NF==0{print $1":"$2"-"$3,$4}' > Brain-specific-peaks.tsv

bedtools intersect -a /tutorial/results/CHIPembryoBrain.narrowPeak -b /tutorial/results/CHIPembryoLiver.narrowPeak -wao|awk 'BEGIN{OFS="\t"; print "Brain_coordinates","Brain_peak"}$NF!=0{print $1,$2,$3,$4}' > Brain-specific-peaks.bed

bedtools intersect -a /tutorial/results/CHIPembryoLiver.narrowPeak -b /tutorial/results/CHIPembryoBrain.narrowPeak -wao|awk 'BEGIN{FS=OFS="\t"; print "Liver_coordinates","Liver_peak"}$NF==0{print $1":"$2"-"$3,$4}' > Liver-specific-peaks.tsv

bedtools intersect -a /tutorial/results/CHIPembryoLiver.narrowPeak -b /tutorial/results/CHIPembryoBrain.narrowPeak -wao|awk 'BEGIN{OFS="\t"; print "Liver_coordinates","Liver_peak"}$NF!=0{print $1,$2,$3,$4}' > Liver-specific-peaks.bed

- BARPLOT -

ls *tsv|while read f; do echo -e $f"\t"$(grep -v nan $f |wc -l);done | ggbarplot.R -i stdin -o number_of_peaks.pdf -f 1 --palette_fill /tutorial/palettes/palTissue.txt --title "H3k4me3 “peaks --y_title "Number of peaks" --x_title "Tissues"

5. Create a BED file of 200bp up/downstream TSS of genes and overlap DEGs (step 1) with the 3 sets of H3K4me3 peaks classified in the previous step (4). Show three
examples in the UCSC genome browser, including RNA-seq, ChIP-seq and ATAC-seq tracks. Ideally, one example of each peak set (i.e. shared peak, peak exclusively called
in brain and peak exclusively called in liver). Discuss the integration of the three datasets in the TSS of the selected cases.

mkdir ../TSS

cd TSS

- subset -

awk '{print $1}' /tutorial/analysis/edgeR.0.01.over_liver_X_brain.txt > /tutorial/TSS/genes_names_liver_X_brain.txt
awk '{print $1}' /tutorial/analysis/edgeR.0.01.over_brain_X_liver.txt > /tutorial/TSS/genes_names_brain_X_liver.txt

grep -Ff genes_names_liver_X_brain.txt protein-coding-genes-200up_downTSS.bed  > TSS_genes_liver_X_brain.bed
grep -Ff genes_names_brain_X_liver.txt protein-coding-genes-200up_downTSS.bed  > TSS_genes_brain_X_liver.bed


- bedtools intersect -

with TSS_genes_liver_X_brain.bed:

1. bedtools intersect -a TSS_genes_liver_X_brain.bed -b /tutorial/chip-analysis/Liver-specific-peaks.bed -wao > TSS-liver_liver_X_brain.bed
2. bedtools intersect -a TSS_genes_liver_X_brain.bed -b /tutorial/chip-analysis/common-peaks-Brain-liver.bed -wao > TSS-liver_common.bed

with  TSS_genes_brain_X_liver.bed:

1. bedtools intersect -a  TSS_genes_brain_X_liver.bed -b /tutorial/chip-analysis/Brain-specific-peaks.bed -wao > TSS-brain_brain_X_liver.bed
2. bedtools intersect -a  TSS_genes_brain_X_liver.bed -b/tutorial/chip-analysis/common-peaks.Brain-liver.bed -wao > TSS-brain_common.bed

