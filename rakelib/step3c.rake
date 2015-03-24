##
## Step 3c - count by coding 5'-UTR & the proximal upstream
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

step3c_bed_sources = Array.new
begin
  conf = CONF[LIBIDS[0]]
  tdbpath = File.expand_path(conf['TRANSCRIPT']+'.txt.gz')
  if /^knownGene/ =~ tdbpath.pathmap('%f')
    step3c_bed_sources.push(tdbpath.sub('knownGene', 'kgXref'))
  else
    step3c_bed.push(tdbPath.sub('ensGene', 'ensemblToGeneName'))
  end
  step3c_bed_sources.push(File.expand_path(conf['GENOMESPIKERIBO']+'.fa.fai'))
  step3c_bed_sources.push(tdbpath)
end

file 'out/cg/regions.bed.gz' => step3c_bed_sources do |t|
  mkdir_p t.name.pathmap('%d')

  acc2sym = Hash.new
  acc2sympath = t.source
  filter = (/^kgXref/ =~ acc2sympath.pathmap('%f')) ? '| gcut -f 1,5' : ''

  infp = open("| gunzip -c #{acc2sympath} #{filter}")
  while line = infp.gets
    acc, sym = line.rstrip.split(/\t/)
    acc2sym[acc] = sym.gsub(' ', '%20')
  end
  infp.close

  ofs = filter != '' ? 1 : 2
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

def step3c_sources(path)
  return ['out/cg/regions.bed.gz',
          path.sub(/^tmp\/cg/, 'tmp').sub(/\.step3c$/, '.step3b')]
end

rule /^tmp\/cg\/[^\/]+\.step3c$/ => [->(path){ step3c_sources(path) }] do |t|
  sym2cnt = Hash.new
  infp = open("| intersectBed -c -s -a #{t.source} -b #{t.sources[1]} | gcut -f 4,7")
  while line = infp.gets
    sym, cnt = line.rstrip.split(/\t/)
    if sym2cnt.key?(sym)
      sym2cnt[sym] += cnt.to_i
    else
      sym2cnt[sym] = cnt.to_i
    end
  end
  infp.close

  mkdir_p t.name.pathmap('%d')
  outfp = open(t.name, 'w')
  sym2cnt.keys.sort.each do |sym|
    outfp.puts [sym, sym2cnt[sym]].join("\t")
  end
  outfp.close
end

task :clean_step3c do
  sh "rm -rf out/cg/regions.bed.gz"
  LIBIDS.each do |libid|
    sh "rm -rf tmp/cg/#{libid}.*.step3c"
  end
end
