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
  #      => [ expression_#{libid}, tmp2 ]
  #
  target = "out/exp/#{libid}.reads.uniq.txt.gz"
  file target => [timestamp]+tmp2 do |t|
    sh "mkdir -p #{t.name.pathmap('%d')}"
    join_counts(t.name, t.prerequisites[1..-1])
  end
  tmp.push(target)
}
task :expression => tmp

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
  sh "bin/_process_overlapping_fixedstep_count.rb #{t.prerequisites.join(' ')} | gzip -c > #{t.name}"
end

####
#
# cleaning
#
task 'clean_expression' do
  LIBIDS.each { |libid| sh "rm -rf tmp/exp/#{libid}.* tmp/#{libid}.expression.timestamp out/exp/#{libid}.*" }
end
