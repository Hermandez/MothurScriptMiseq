#!/bin/bash
cd /home/Git/armel/Armel16SsequencesBackup/working/
#-- Set working directory
mothur '#make.file(inputdir=., type=gz, prefix=armel16s)'

# Reducing sequencing and PCR errors
mothur '#make.contigs(file=armel16s.files, processors=8)' # Make armel16s.files containing group names and paired reads (NB: Make sure your path contains no underscore)
mothur '#summary.seqs(fasta=armel16s.trim.contigs.fasta)' # Examine your read lengths and pick an appropraite maxlength for the next step
mothur '#screen.seqs(fasta=armel16s.trim.contigs.fasta, group=armel16s.contigs.groups, maxambig=0, maxlength=480)'
#mothur '#screen.seqs(fasta=armel16s.trim.contigs.fasta, group=armel16s.contigs.groups, summary=armel16s.trim.contigs.summary, maxambig=0, maxlength=470)'
#mothur '#get.current()'
mothur '#summary.seqs(fasta=armel16s.trim.contigs.good.fasta)'
#mothur '#summary.seqs(fasta=current)'
#mothur '#summary.seqs()'

# Processing improved sequences
mothur '#unique.seqs(fasta=armel16s.trim.contigs.good.fasta)' # Assuming there are duplicates reads, take only uniq
mothur '#count.seqs(name=armel16s.trim.contigs.good.names, group=armel16s.contigs.good.groups)'
mothur '#summary.seqs(count=armel16s.trim.contigs.good.count_table)'
mothur '#pcr.seqs(fasta=silva.bacteria.fasta, start=5650, end=25519, keepdots=F, processors=8)' # Set the right coordinates to cover your V3 V4 region
mothur '#rename.file(input=silva.bacteria.pcr.fasta, new=silva.v4.fasta)' # Rename the truncated reference
mothur '#summary.seqs(fasta=silva.v4.fasta)'
mothur '#align.seqs(fasta=armel16s.trim.contigs.good.unique.fasta, reference=silva.v4.fasta)' # Now aling your reads to the reference
mothur '#summary.seqs(fasta=armel16s.trim.contigs.good.unique.align, count=armel16s.trim.contigs.good.count_table)'
mothur '#screen.seqs(fasta=armel16s.trim.contigs.good.unique.align, count=armel16s.trim.contigs.good.count_table, summary=armel16s.trim.contigs.good.unique.summary, start=738, end=19666, maxhomop=8)' # Check the read leangths and their start and end positions from the summary above. 
mothur '#summary.seqs(fasta=current, count=current)'
mothur '#filter.seqs(fasta=armel16s.trim.contigs.good.unique.good.align, vertical=T, trump=.)'
mothur '#unique.seqs(fasta=armel16s.trim.contigs.good.unique.good.filter.fasta, count=armel16s.trim.contigs.good.good.count_table)'
mothur '#pre.cluster(fasta=armel16s.trim.contigs.good.unique.good.filter.unique.fasta, count=armel16s.trim.contigs.good.unique.good.filter.count_table, diffs=2)'
mothur '#chimera.vsearch(fasta=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)' # Search for chimeras. Notice the chemras found with and without taking into account abundance
mothur '#remove.seqs(fasta=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)'
mothur '#summary.seqs(fasta=current, count=current)'
mothur '#classify.seqs(fasta=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=trainset9_032012.pds.fasta, taxonomy=trainset9_032012.pds.tax, cutoff=80)'
mothur '#remove.lineage(fasta=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)'
mothur '#summary.tax(taxonomy=current, count=current)'

# Assessing error rates (Run this only if you co-sequnced your samples with the Mock community of any other reference)
#mothur '#get.groups(count=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, fasta=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, groups=Mock)'
#mothur '#seq.error(fasta=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta, count=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count_table, reference=HMP_MOCK.v35.fasta, aligned=F)'

# Cluster sequences into OTUs and check spurious OTUs
mothur '#dist.seqs(fasta=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, cutoff=0.03)'
mothur '#cluster(column=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.dist, count=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, cutoff=0.03)'
mothur '#make.shared(list=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.list, count=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, label=0.03)'
mothur '#rarefaction.single(shared=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.shared)'

# Preparing for analysis
#mothur '#remove.groups(count=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, fasta=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, taxonomy=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, groups=Mock)'

# OTUs
mothur '#dist.seqs(fasta=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, cutoff=0.03)'
mothur '#cluster(column=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.dist, count=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, cutoff=0.03)'
mothur '#cluster.split(fasta=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, splitmethod=classify, taxlevel=4, cutoff=0.03)'
mothur '#make.shared(list=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.list, count=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, label=0.03)'
mothur '#classify.otu(list=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.list, count=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, label=0.03)'

# Phylotypes
mothur '#phylotype(taxonomy=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy)'
mothur '#make.shared(list=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.list, count=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, label=1)'
mothur '#classify.otu(list=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.list, count=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, label=1)'

# Phylogenetic
mothur '#dist.seqs(fasta=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, output=lt, processors=8)'
mothur '#clearcut(phylip=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.phylip.dist)'

# Analysis
mothur '#rename.file(taxonomy=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.0.03.cons.taxonomy, shared=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.shared)'
mothur '#count.groups(shared=armel16s.opti_mcc.shared)' # Check the sequence size of the samples (groups)
mothur '#sub.sample(shared=armel16s.opti_mcc.shared, size=6347)' # Use the smallest sequence size 

# OTU-based analysis
mothur '#rarefaction.single(shared=armel16s.opti_mcc.shared, calc=sobs, freq=100)'
mothur '#summary.single(shared=armel16s.opti_mcc.shared, calc=nseqs-coverage-sobs-invsimpson, subsample=T)'

# Beta diversity measurements
mothur '#dist.shared(shared=armel16s.opti_mcc.shared, calc=thetayc-jclass, subsample=t)'
mothur '#pcoa(phylip=armel16s.opti_mcc.thetayc.0.03.lt.ave.dist)'
mothur '#nmds(phylip=armel16s.opti_mcc.thetayc.0.03.lt.ave.dist)'
mothur '#nmds(phylip=armel16s.opti_mcc.thetayc.0.03.lt.ave.dist, mindim=3, maxdim=3)'
#mothur '#amova(phylip=armel16s.opti_mcc.thetayc.0.03.lt.ave.dist, design=mouse.time.design)' # You'll need your metadata here to create a file with groups (samples) in first column and the variable to measure in the second column with no header. Replace mouse.time.design with the new file
#mothur '#homova(phylip=armel16s.opti_mcc.thetayc.0.03.lt.ave.dist, design=mouse.time.design)' # Same thing as above. You can therefore run multiple amovas with as many variables you have
mothur '#corr.axes(axes=armel16s.opti_mcc.thetayc.0.03.lt.ave.pcoa.axes, shared=armel16s.opti_mcc.0.03.subsample.shared, method=spearman, numaxes=3)'
#mothur '#corr.axes(axes=armel16s.opti_mcc.thetayc.0.03.lt.ave.pcoa.axes, metadata=mouse.dpw.metadata, method=spearman, numaxes=3)' # This requires that you provide data indicating metadata about each sample. For example, you may know the weight, height, blood pressure, etc. of the subjects in the samples. The file has two columns, col1 is the group (sample) name and col2 the variable with a header as: group variable
mothur '#get.communitytype(shared=armel16s.opti_mcc.0.03.subsample.shared)'

# Population-level analysis (Again, you need the metadata with variables to run these lines)
#mothur '#metastats(shared=armel16s.opti_mcc.0.03.subsample.shared, design=mouse.time.design)' #groups=F3D0-F3D1-F3D141-F3D142-F3D144-F3D145-F3D146-F3D147-F3D148-F3D149-F3D150-F3D2-F3D3-F3D5-F3D6-F3D7-F3D8-F3D9, 
#mothur '#lefse(shared=armel16s.opti_mcc.0.03.subsample.shared, design=mouse.time.design)'

# Phylogeny-based analysis
# Alpha diversity
mothur '#phylo.diversity(tree=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.phylip.tre, count=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, rarefy=T)'
mothur '#unifrac.unweighted(tree=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.phylip.tre, count=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, distance=lt, processors=2, random=F, subsample=t)'
mothur '#unifrac.weighted(tree=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.phylip.tre, count=armel16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, distance=lt, processors=2, random=F, subsample=t)'

