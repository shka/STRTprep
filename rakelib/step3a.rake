##
## Step 3a - Region of coding 5'-UTR & the proximal upstream
##

def extract_5utr(acc2sym, outfp, tbl, chrEnds, ofs=1)
  infp = open("| gunzip -c #{tbl} | gcut -f #{ofs}-")
  while line = infp.gets
    cols = line.rstrip.split(/\t/)
    # no ORF ~ CDS start != CDS stop, since the start position is 0-based
    if cols[5] != cols[6] && chrEnds.key?(cols[1])
      lefts = cols[8].split(/,/)
      rights = cols[9].split(/,/)
      sym = acc2sym.key?(cols[0]) ? acc2sym[cols[0]] : cols[11]
      if cols[2] == '+'
        cdsleft = cols[5].to_i
        0.upto(lefts.length-1) do |i|
          left = lefts[i].to_i
          right = rights[i].to_i
          if right < cdsleft
            outfp.puts [cols[1], left, right, sym, 0, '+'].join("\t")
          elsif left < cdsleft
            outfp.puts [cols[1], left, cdsleft, sym, 0, '+'].join("\t")
          end
        end
      else
        cdsright = cols[6].to_i
        (rights.length-1).downto(0) do |i|
          left = lefts[i].to_i
          right = rights[i].to_i
          if cdsright < left
            outfp.puts [cols[1], left, right, sym, 0, '-'].join("\t")
          elsif cdsright < right
            outfp.puts [cols[1], cdsright, right, sym, 0, '-'].join("\t")
          end
        end
      end
    end
  end
  infp.close
end

def extract_proxup(acc2sym, outfp, tbl, chrEnds, ofs=1)
  infp = open("| gunzip -c #{tbl} | gcut -f #{ofs}-")
  while line = infp.gets
    cols = line.rstrip.split(/\t/)
    # no ORF ~ CDS start != CDS stop, since the start position is 0-based
    if cols[5] != cols[6] && chrEnds.key?(cols[1])
      lefts = cols[8].split(/,/)
      rights = cols[9].split(/,/)
      sym = acc2sym.key?(cols[0]) ? acc2sym[cols[0]] : cols[11]
      if cols[2] == '+'
        left = lefts[0].to_i
        outfp.puts [cols[1], (left-500 > 0 ? left-500: 0), left, sym, 0, '+'].join("\t")
      else
        right = rights[-1].to_i
        outfp.puts [cols[1], right, (right+500 > chrEnds[cols[1]] ? chrEnds[cols[1]] : right+500), sym, 0, '-'].join("\t")
      end
    end
  end
  infp.close
end

step3a_bed_sources = Array.new
begin
  tdbpath = File.expand_path(PREPROCESS['TRANSCRIPT'].pathmap('%d'))
  ref = File.expand_path(PREPROCESS['GENOMESPIKERIBO']+'.fa.fai')
  if File.exist?("#{tdbpath}/kgXref.txt.gz")
    step3a_bed_sources.push("#{tdbpath}/kgXref.txt.gz")
    step3a_bed_sources.push(ref)
    step3a_bed_sources.push("#{tdbpath}/knownGene.txt.gz")
  elsif File.exist?("#{tdbpath}/refGene.txt.gz")
    step3a_bed_sources.push("#{tdbpath}/refGene.txt.gz")
    step3a_bed_sources.push(ref)
    step3a_bed_sources.push("#{tdbpath}/refGene.txt.gz")
  elsif File.exist?("#{tdbpath}/ensemblToGeneName.txt.gz")
    step3a_bed_sources.push("#{tdbpath}/ensemblToGeneName.txt.gz")
    step3a_bed_sources.push(ref)
    step3a_bed_sources.push("#{tdbpath}/ensGene.txt.gz")
  else
    raise "No annotation file at #{tdbpath}."
  end
end

def load_acc2sym(path)
  acc2sym = Hash.new
  cut = ''
  cut = '| gcut -f 1,5' if /^kgXref/ =~ path.pathmap('%f')
  cut = '| gcut -f 2,13' if /^refGene/ =~ path.pathmap('%f')
  infp = open("| gunzip -c #{path} #{cut}")
  while line = infp.gets
    acc, sym = line.rstrip.split(/\t/)
    acc2sym[acc] = sym.gsub(' ', '%20')
  end
  infp.close

  return acc2sym
end

file 'out/byGene/regions.bed.gz' => step3a_bed_sources do |t|
  mkdir_p t.name.pathmap('%d')

  acc2sym = load_acc2sym(t.source)

  outfp = open("| gsort -t '\t' -k 1,1 -k 2,2n | mergeBed -s -o distinct,distinct,distinct -c 4,5,6 -i - | grep -v , | gzip -c > #{t.name}", 'w')

  chrEnds = Hash.new
  infp = open(t.sources[1])
  while line = infp.gets
    cols = line.rstrip.split(/\t/)
    outfp.puts [cols[0], 0, 50, cols[0], 0, '+'].join("\t") if /^RNA_SPIKE_/ =~ cols[0]
    chrEnds[cols[0]] = cols[1].to_i
  end
  infp.close

  tbl = t.sources[t.sources.length == 2 ? 0 : 2]
  ofs = (/^kgXref/ =~ t.source.pathmap('%f')) ? 1 : 2
  extract_5utr(acc2sym, outfp, tbl, chrEnds, ofs)
  extract_proxup(acc2sym, outfp, tbl, chrEnds, ofs)
  outfp.close
end

#

def extract_exon(acc2sym, outfp, tbl, chrEnds, ofs=1)
  infp = open("| gunzip -c #{tbl} | gcut -f #{ofs}-")
  while line = infp.gets
    cols = line.rstrip.split(/\t/)
    if cols[5] != cols[6] && chrEnds.key?(cols[1])
      lefts = cols[8].split(/,/)
      rights = cols[9].split(/,/)
      0.upto(lefts.length-1) do |i|
        left = lefts[i].to_i
        right = rights[i].to_i
        sym = acc2sym.key?(cols[0]) ? acc2sym[cols[0]] : cols[11]
        outfp.puts [cols[1], left, right, sym, 0, cols[2]].join("\t")
      end
    end
  end
  infp.close
end

file 'tmp/byGene/regions_forQC.bed.gz' => step3a_bed_sources do |t|
  mkdir_p t.name.pathmap('%d')

  acc2sym = load_acc2sym(t.source)
  
  outfp = open("| gsort -t '\t' -k 1,1 -k 2,2n | mergeBed -s -o distinct,distinct,distinct -c 4,5,6 -i - | gzip -c > #{t.name}", 'w')

  chrEnds = Hash.new
  infp = open(t.sources[1])
  while line = infp.gets
    cols = line.rstrip.split(/\t/)
    outfp.puts [cols[0], 0, 50, cols[0], 0, '+'].join("\t") if /^RNA_SPIKE_/ =~ cols[0]
    chrEnds[cols[0]] = cols[1].to_i
  end
  infp.close

  tbl = t.sources[t.sources.length == 2 ? 0 : 2]
  ofs = (/^kgXref/ =~ t.source.pathmap('%f')) ? 1 : 2
  extract_exon(acc2sym, outfp, tbl, chrEnds, ofs)
  extract_proxup(acc2sym, outfp, tbl, chrEnds, ofs)
  outfp.close
end
