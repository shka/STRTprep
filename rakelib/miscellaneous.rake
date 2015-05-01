require 'mkfifo'
require 'tempfile'

@mytemppaths = Array.new

def mytemppath(basename, tmpdir = Dir::tmpdir)
  fp = Tempfile.open(basename, tmpdir)
  path = fp.path
  @mytemppaths.push(path)
  fp.close!
  path
end

END { @mytemppaths.each { |path| sh "rm -rf #{path}" if File.exist?(path) } }

def mymkfifo(basename, tmpdir = Dir::tmpdir)
  path = mytemppath(basename, tmpdir)
  File.mkfifo(path)
  path
end

def parse_librunid(path)
  librunid = path.pathmap('%n')
  return librunid.pathmap('%n'), librunid.pathmap('%x')[1..-1].to_i
end
