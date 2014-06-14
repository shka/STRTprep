require 'parallel'

####
#
# dummy paths
#
tmp = Array.new
LIBIDS.each { |libid|
  tmp2 = Array.new
  open(File.expand_path(CONF[libid]['LAYOUT'])).each { |line|
    cols = line.rstrip.split
    tmp2.push("tmp/exp/#{libid}.#{cols[0]}.txt.gz")
  }
  #
  timestamp = "tmp/#{libid}.expression.timestamp"
  file timestamp => tmp2 do |t|
    sh "touch #{t.name}"
    begin
      Parallel.map(t.prerequisites, :in_threads => PROCS) { |target|
        Rake::Task[target].invoke
      }
    rescue
      sh "rm -rf #{timestamp}"
    end
  end
  tmp.push(timestamp)
  ####
  #
  # file out/exp/#{libid}.reads.uniq.txt.gz
  #      => [ tmp/#{libid}.expression.timestamp, tmp2 ]
  #
  table = "out/exp/#{libid}.reads.uniq.txt.gz"
  file table => [timestamp]+tmp2 do |t|
    sh "mkdir -p out/exp"
    join_counts(t.name, t.prerequisites[1..-1])
  end
  #
  tmp.push("out/exp/#{libid}.pvclust.uniq.success.fluctuated.RData.gz",
           "out/exp/#{libid}.pca.uniq.success.fluctuated.RData.gz",
           "out/exp/#{libid}.annotation.uniq.txt.gz")
}
task :expression => tmp + [:hub]

def join_counts(out, ins)
  reg2cnts = Hash.new
  colnames = Array.new
  cols = ins.length
  ins_sorted = ins.sort
  ins.each { |target|
    i = ins_sorted.index(target)
    colnames.push(/exp\/([^.]+\.[^\.]+)\.txt/.match(target).to_a[1])
    begin
      open("| gunzip -c #{target}").each { |line|
        reg, cnt = line.rstrip.split(/\t/)
        reg2cnts[reg] = Array.new(cols) { 0 } unless reg2cnts.key?(reg)
        reg2cnts[reg][i] = cnt
      }
    rescue
      puts "... no alignment in #{target}; skipped"
    end
  }
  fp = open("| gzip -c --best > #{out}", 'w')
  fp.puts (['']+colnames).join("\t")
  reg2cnts.keys.sort.each { |reg| fp.puts ([reg]+reg2cnts[reg]).join("\t") }
  fp.close
end

####
#
# file tmp/exp/#{libid}.bed.gz
#      => CONF[#{libid}]['GENOMESPIKERIBO']+'.fa.fai'
#
rule /tmp\/exp\/[^.]+\.bed.gz/ => proc { |target|
  File.expand_path(CONF[/exp\/([^.]+)\./.match(target).to_a[1]]['GENOMESPIKERIBO']+'.fa.fai')
} do |t|
  sh "mkdir -p #{t.name.pathmap('%d')}"
  fp = open("| gzip -c > #{t.name}", 'w')
  open("| cut -f 1,2 #{t.prerequisites[0]}").each { |line|
    chr, size = line.rstrip.split(/\t/)
    pos = 1
    while pos < size.to_i
      fp.puts [chr, pos-1, pos+49, "#{chr}:#{pos}-#{pos+49},+", '.', '+'].join("\t")
      fp.puts [chr, pos-1, pos+49, "#{chr}:#{pos}-#{pos+49},-", '.', '-'].join("\t")
      pos += 25
    end
  }
  fp.close
end

####
#
# file tmp/exp/#{libid}.#{wellid}.txt.gz
#      => [ tmp/hub/#{libid}.#{wellid}.fwd.uniq.5p.bed.gz,
#           tmp/hub/#{libid}.#{wellid}.rev.uniq.5p.bed.gz ]
#
rule /tmp\/exp\/[^.]+\.[^.]+\.txt/ => proc { |target|
  [ target.sub('exp', 'hub').sub('.txt.gz', '.fwd.uniq.5p.bed.gz'),
    target.sub('exp', 'hub').sub('.txt.gz', '.rev.uniq.5p.bed.gz') ]
} do |t|
  sh "mkdir -p tmp/exp"
  sh "bin/_process_overlapping_fixedstep_count.rb #{t.prerequisites.join(' ')} | gzip -c > #{t.name}"
end

####
#
# file out/exp/#{libid}.pvclust.uniq.success.fluctuated.RData.gz
#      => out/exp/#{libid}.nreads.uniq.success.fluctuated.RData.gz
#
rule /\.pvclust\.uniq\.success\.fluctuated\.RData\.gz/ => proc { |target|
  target.sub('.pvclust.', '.nreads.')
} do |t|
  sh "R --vanilla --quiet --args #{t.prerequisites[0]} #{PROCS} < bin/_process_expression_pvclust.R"
end

####
#
# file out/exp/#{libid}.pca.uniq.success.fluctuated.RData.gz
#      => out/exp/#{libid}.nreads.uniq.success.fluctuated.RData.gz
#
rule /\.pca\.uniq\.success\.fluctuated\.RData\.gz/ => proc { |target|
  target.sub('.pca.', '.nreads.')
} do |t|
  sh "R --vanilla --quiet --args #{t.prerequisites[0]} < bin/_process_expression_pca.R"
end

####
#
# file out/exp/#{libid}.nreads.uniq.success.fluctuated.RData.gz
#      => out/exp/#{libid}.nreads.uniq.success.RData.gz
#
rule /\.nreads\.uniq\.success\.fluctuated\.RData\.gz/ => proc { |target|
  target.sub('.fluctuated', '')
} do |t|
  libid = /\/([^\/.]+)\.nreads/.match(t.name).to_a[1]
  sh "R --vanilla --quiet --args #{t.prerequisites[0]} '#{CONF[libid]['WELLS']['EXCEPTION']}' #{CONF[libid]['FDR']} < bin/_process_expression_fluctuation.R"
end

####
#
# file out/exp/#{libid}.nreads.uniq.success.RData.gz
#      => out/exp/#{libid}.reads.uniq.txt.gz
#
rule /\.nreads\.uniq\.success\.RData\.gz/ => proc { |target|
  target.sub('.success.RData', '.txt').sub('.nreads', '.reads')
} do |t|
  libid = /\/([^\/.]+)\.nreads/.match(t.name).to_a[1]
  sh "R --vanilla --quiet --args #{t.prerequisites[0]} '#{CONF[libid]['WELLS']['FAILURE']}' '#{CONF[libid]['WELLS']['EXCEPTION']}' < bin/_process_expression.R"
end

####
#
# file tmp/exp/#{libid}.reads.uniq.bed.gz
#      => out/exp/#{libid}.reads.uniq.txt.gz
#
rule /[^.]+\.reads\.uniq\.bed\.gz/ => proc { |target|
  target.sub('.bed', '.txt').sub('tmp/', 'out/')
} do |t|
  sh "gunzip -c #{t.prerequisites[0]} | cut -f 1 | grep -v -P '^$' | ruby -nle 't=/([^:]+):(\\d+)\-(\\d+),([-+])/.match($_).to_a; puts [t[1], t[2].to_i-1, t[3], t[0], 0, t[4]].join(\"\\t\")' | sort -k1,1 -k2,2n | gzip -c > #{t.name}"
end

####
#
# file tmp/exp/#{libid}.reads.uniq.class\d.txt.gz
#      => [ tmp/exp/#{libid}.reads.uniq.bed.gz, tmp/*.class\d.bed.gz ]
#
rule /[^.]+\.reads\.uniq\.class\d\.txt\.gz/ => proc { |target|
  libid = /\/([^\/.]+)\.reads/.match(target).to_a[1]
  classid = /\.class(\d)/.match(target).to_a[1]
  [ target.sub("\.class#{classid}.txt", '.bed'),
    File.expand_path(CONF[libid]['TRANSCRIPT'] + ".class#{classid}.bed.gz") ]
} do |t|
  sh "bedtools intersect -wa -wb -s -a #{t.prerequisites[0]} -b #{t.prerequisites[1]} | cut -f 4,10 | cut -d '|' -f 1 | sort -u | gzip -c > #{t.name}"
end

###
#
# file out/exp/#{libid}.annotation.uniq.txt.gz
#      => [ out/exp/#{libid}.reads.uniq.txt.gz,
#           tmp/exp/#{libid}.reads.uniq.class\d.txt.gz ]
#
rule /[^.]+\.annotation\.uniq\.txt\.gz/ => proc {  |target|
  tmp = [ target.sub('.annotation', '.reads') ]
  1.upto(8) { |i| tmp.push(target.sub('out/', 'tmp/').sub('.annotation.uniq', ".reads.uniq.class#{i}")) }
  tmp
} do |t|
  anncls = ['Unannotated', "Coding 5'-UTR", "Coding upstream", 'Coding CDS', "Coding 3'-UTR", 'Noncoding 1st-exon', 'Noncoding upstream', 'Noncoding other exon', 'Intron']
  reg2cls = Hash.new
  reg2sym = Hash.new
  1.upto(8) { |i|
    open("| gunzip -c #{t.prerequisites[i]} ").each { |line|
      reg, sym = line.rstrip.split(/\t/)
      reg2cls[reg] = anncls[i] unless reg2cls.key?(reg)
      if reg2cls[reg] == anncls[i]
        reg2sym[reg] = Hash.new unless reg2sym.key?(reg)
        reg2sym[reg][sym] = ''
      end
    }
  }
  fp = open("| gzip -c > #{t.name}", 'w')
  open("| gunzip -c #{t.prerequisites[0]} | cut -f 1").each { |line|
    reg = line.rstrip
    if reg == ''
      fp.puts ['', 'Symbol', 'Location'].join("\t")
    else
      cls = reg2cls.key?(reg) ? reg2cls[reg] : 'Unannotated'
      sym = reg2cls.key?(reg) ? reg2sym[reg].keys.join(',') : reg
      fp.puts [reg, sym, cls].join("\t")
    end
  }
  fp.close
end

####
#
# file out/exp/#{libid}.reads.uniq.detected.normalized.fluctuated.annotated.txt.gz
#      => [ out/exp/#{libid}.reads.uniq.detected.normalized.fluctuated.RData.gz,
#           out/exp/#{libid}.annotation.uniq.txt.gz ]
#
rule /\.reads\.uniq\.detected\.normalized\.fluctuated\.annotated\.txt\.gz$/ => proc { |target|
  [ target.sub('.annotated.txt', '.RData'),
    target.sub('.reads', '.annotation').sub('.detected.normalized.fluctuated.annotated', '')] } do |t|
  sh "R --vanilla --quiet --args #{t.prerequisites.join("\t")} #{t.name} < bin/_process_annotation.R"
end

####
#
# cleaning
#
task 'clean_expression' do
  LIBIDS.each { |libid| sh "rm -rf tmp/exp/#{libid}.* tmp/#{libid}.expression.timestamp out/exp/#{libid}.*" }
end
