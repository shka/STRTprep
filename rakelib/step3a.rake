##
## Step 3a - Region of coding 5'-UTR & the proximal upstream
##

def extract_5utr(acc2sym, outfp, tbl, ofs=1)
  infp = open("| gunzip -c #{tbl} | cut -f #{ofs}-")
  while line = infp.gets
    cols = line.rstrip.split(/\t/)
    # no ORF ~ CDS start != CDS stop, since the start position is 0-based
    if cols[5] != cols[6]
      lefts = cols[8].split(/,/)
      rights = cols[9].split(/,/)
      if cols[2] == '+'
        cdsleft = cols[5].to_i
        0.upto(lefts.length-1) do |i|
          left = lefts[i].to_i
          right = rights[i].to_i
          if right < cdsleft
            outfp.puts [cols[1], left, right, acc2sym[cols[0]], 0, '+'].join("\t")
          elsif left < cdsleft
            outfp.puts [cols[1], left, cdsleft, acc2sym[cols[0]], 0, '+'].join("\t")
          end
        end
      else
        cdsright = cols[6].to_i
        (rights.length-1).downto(0) do |i|
          left = lefts[i].to_i
          right = rights[i].to_i
          if cdsright < left
            outfp.puts [cols[1], left, right, acc2sym[cols[0]], 0, '-'].join("\t")
          elsif cdsright < right
            outfp.puts [cols[1], cdsright, right, acc2sym[cols[0]], 0, '-'].join("\t")
          end
        end
      end
    end
  end
  infp.close
end

def extract_proxup(acc2sym, outfp, tbl, ofs=1)
  infp = open("| gunzip -c #{tbl} | cut -f #{ofs}-")
  while line = infp.gets
    cols = line.rstrip.split(/\t/)
    # no ORF ~ CDS start != CDS stop, since the start position is 0-based
    if cols[5] != cols[6]
      lefts = cols[8].split(/,/)
      rights = cols[9].split(/,/)
      if cols[2] == '+'
        left = lefts[0].to_i
        outfp.puts [cols[1], left-500, left, acc2sym[cols[0]], 0, '+'].join("\t")
      else
        right = rights[-1].to_i
        outfp.puts [cols[1], right, right+500, acc2sym[cols[0]], 0, '-'].join("\t")
      end
    end
  end
  infp.close
end

step3a_bed_sources = Array.new
begin
  conf = CONF[LIBIDS[0]]
  tdbpath = File.expand_path(conf['TRANSCRIPT']+'.txt.gz')
  if /^knownGene/ =~ tdbpath.pathmap('%f')
    step3a_bed_sources.push(tdbpath.sub('knownGene', 'kgXref'))
  else
    step3a_bed_sources.push(tdbpath.sub('ensGene', 'ensemblToGeneName'))
  end
  step3a_bed_sources.push(File.expand_path(conf['GENOMESPIKERIBO']+'.fa.fai'))
  step3a_bed_sources.push(tdbpath)
end

def load_acc2sym(path, isKnownGene)
  acc2sym = Hash.new
  infp = open("| gunzip -c #{path} #{isKnownGene ? '| gcut -f 1,5' : ''}")
  while line = infp.gets
    acc, sym = line.rstrip.split(/\t/)
    acc2sym[acc] = sym.gsub(' ', '%20')
  end
  infp.close

  return acc2sym
end

file 'out/cg/regions.bed.gz' => step3a_bed_sources do |t|
  mkdir_p t.name.pathmap('%d')

  isKnownGene = /^kgXref/ =~ t.source.pathmap('%f')

  acc2sym = load_acc2sym(t.source, isKnownGene)
  
  ofs = isKnownGene ? 1 : 2
  outfp = open("| gsort -k 1,1 -k 2,2n | mergeBed -s -o distinct,distinct,distinct -c 4,5,6 -i - | grep -v , | gzip -c > #{t.name}", 'w')
    
  infp = open("| grep ^RNA_SPIKE_ #{t.sources[1]} | cut -f 1")
  while line = infp.gets
    outfp.puts [line.rstrip, 0, 50, line.rstrip, 0, '+'].join("\t")
  end
  infp.close
  
  extract_5utr(acc2sym, outfp, t.sources[2], ofs)
  extract_proxup(acc2sym, outfp, t.sources[2], ofs)
  outfp.close
end

#

def extract_exon(acc2sym, outfp, tbl, ofs=1)
  infp = open("| gunzip -c #{tbl} | cut -f #{ofs}-")
  while line = infp.gets
    cols = line.rstrip.split(/\t/)
    if cols[5] != cols[6]
      lefts = cols[8].split(/,/)
      rights = cols[9].split(/,/)
      0.upto(lefts.length-1) do |i|
        left = lefts[i].to_i
        right = rights[i].to_i
        outfp.puts [cols[1], left, right, acc2sym[cols[0]], 0, cols[2]].join("\t")
      end
    end
  end
  infp.close
end

file 'tmp/cg/regions_forQC.bed.gz' => step3a_bed_sources do |t|
  mkdir_p t.name.pathmap('%d')

  isKnownGene = /^kgXref/ =~ t.source.pathmap('%f')

  acc2sym = load_acc2sym(t.source, isKnownGene)
  
  ofs = isKnownGene ? 1 : 2
  outfp = open("| gsort -k 1,1 -k 2,2n | mergeBed -s -o distinct,distinct,distinct -c 4,5,6 -i - | gzip -c > #{t.name}", 'w')

  infp = open("| grep ^RNA_SPIKE_ #{t.sources[1]} | cut -f 1")
  while line = infp.gets
    outfp.puts [line.rstrip, 0, 50, line.rstrip, 0, '+'].join("\t")
  end
  infp.close
  
  extract_exon(acc2sym, outfp, t.sources[2], ofs)
  extract_proxup(acc2sym, outfp, t.sources[2], ofs)
  outfp.close
end

#

task :clean_step3a do
  sh 'rm out/cg/regions.bed.gz'
  sh 'rm tmp/cg/regions_forQC.bed.gz'
end
