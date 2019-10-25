#bowtie2 --very-sensitive-local -p 2 -x hgr_5s -1 na19240_pcr_10m_1.fq.gz -2 na19240_pcr_10m_2.fq.gz | samtools view -bS - > hgr_5s.test.bam

#bowtie2 --very-sensitive-local -p 2 -x hgr0_draft -1 na19240_pcr_10m_1.fq.gz -2 na19240_pcr_10m_2.fq.gz | samtools view -bS - > hgr0_draft.test.bam

#bowtie2 --very-sensitive-local -p 2 -x 45S_5S_ref -1 na19240_pcr_10m_1.fq.gz -2 na19240_pcr_10m_2.fq.gz | samtools view -bS - > 45S_5S_ref.test.bam

#bowtie2 --very-sensitive-local -p 2 -x hgr0_draft2 -1 na19240_pcr_10m_1.fq.gz -2 na19240_pcr_10m_2.fq.gz | samtools view -bS - > hgr0_draft2.test.bam

#bowtie2 --very-sensitive-local -p 2 -x hgr0_draft3 -1 na19240_pcr_10m_1.fq.gz -2 na19240_pcr_10m_2.fq.gz | samtools view -bS - > hgr0_draft3.test.bam

#bowtie2 --very-sensitive-local -p 2 -x hgr0_draft4 -1 na19240_pcr_10m_1.fq.gz -2 na19240_pcr_10m_2.fq.gz | samtools view -bS - > hgr0_draft4.test.bam
