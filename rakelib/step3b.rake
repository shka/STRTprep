##
## Step 3b - Count in the regions
##

def step3b_sources(path)
  return ['out/byGene/regions.bed.gz',
          path.sub(/^tmp\/byGene/, 'tmp').sub(/\.step3b$/, '.step2e')]
end

def step3b_job(t)
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

rule /^tmp\/byGene\/[^\/]+\.step3b$/ => [->(path){ step3b_sources(path) }] do |t|
  step3b_job(t)
end

#

def step3b_cnt_job(t)
  cnt_coding = 0
  cnt_spike = 0
  infp = open(t.source)
  while line = infp.gets
    acc, cnt = line.rstrip.split(/\t/)
    if /^RNA_SPIKE_/ =~ acc
      cnt_spike += cnt.to_i
    else
      cnt_coding += cnt.to_i
    end
  end
  infp.close
  
  outfp = open(t.name, 'w')
  outfp.puts [cnt_coding, cnt_spike].join("\t")
  outfp.close
end

rule '.step3b_cnt' => '.step3b' do |t|
  step3b_cnt_job(t)
end
