require 'rake'
require 'fileutils'

CSBLAST_VERSION = ENV['CSBLAST_VERSION'] || %x{egrep 'kVersionNumber =' src/application.cc | sed -r 's/^.+\"([0-9]\.[0-9]\.[0-9])\";$/\\1/'}.strip
RELEASE_DIR = "CS-BLAST-" + CSBLAST_VERSION
MAC         = "macmini"
WEBSERVER   = "ws01"
FTP_DIR     = "/home/ftp/data/csblast"
FILES       = FileList['CHANGELOG', 'README', 'LICENSE', 'data/K4000.lib']

task :default => :all

desc "Build tar archives for all architectures"
task :all => [:linux32, :linux64, :macosx]

desc "Remove generated files"
task :clean do
  sh "cd src; make clean"
  rm_rf(RELEASE_DIR)
  sh "ssh #{MAC} 'rm -rf src/cs'"
end

desc "Build tar archive for Linux 32bit"
task :linux32 do
  dirname = File.join(RELEASE_DIR, "csblast-#{CSBLAST_VERSION}-linux32")
  setup_target_dir dirname

  sh "cd src; make clean; make csblast csbuild_aa NDEBUG=1 M32=1"
  mv "src/csblast", dirname
  mv "src/csbuild_aa", File.join(dirname, "csbuild")

  sh "cd src; make clean; make csblast csbuild_aa NDEBUG=1 M32=1 STATIC=1"
  mv "src/csblast", File.join(dirname, "csblast_static")
  mv "src/csbuild_aa", File.join(dirname, "csbuild_static")

  tar_and_copy_to_webserver dirname
end

desc "Build tar archive for Linux 64bit"
task :linux64 do
  dirname = File.join(RELEASE_DIR, "csblast-#{CSBLAST_VERSION}-linux64")
  setup_target_dir dirname

  sh "cd src; make clean; make csblast csbuild_aa NDEBUG=1"
  mv "src/csblast", dirname
  mv "src/csbuild_aa", File.join(dirname, "csbuild")

  sh "cd src; make clean; make csblast csbuild_aa NDEBUG=1 STATIC=1"
  mv "src/csblast", File.join(dirname, "csblast_static")
  mv "src/csbuild_aa", File.join(dirname, "csbuild_static")

  tar_and_copy_to_webserver dirname
end

desc "Build tar archive for Mac OS X"
task :macosx do
  dirname = File.join(RELEASE_DIR, "csblast-#{CSBLAST_VERSION}-macosx")
  setup_target_dir dirname

  sh "rsync -avzrq --delete --exclude=.git --exclude=build --exclude=data --exclude=#{RELEASE_DIR} ./ #{MAC}:src/cs"

  sh "ssh #{MAC} 'cd src/cs/src; make clean; make csblast csbuild_aa NDEBUG=1'"
  sh "scp #{MAC}:src/cs/src/csblast #{dirname}/csblast"
  sh "scp #{MAC}:src/cs/src/csbuild_aa #{dirname}/csbuild"

  sh "ssh #{MAC} 'rm -rf src/cs'"

  tar_and_copy_to_webserver dirname
end

def setup_target_dir dirname
  mkdir_p dirname
  FILES.each do |f|
    cp f, dirname
  end
end

def tar_and_copy_to_webserver dirname
  basedir  = File.dirname dirname
  basename = File.basename dirname

  sh "cd #{basedir}; tar czf #{basename}.tar.gz #{basename}"
  sh "ssh #{WEBSERVER} 'mkdir -p #{FTP_DIR}/#{RELEASE_DIR}'"
  sh "scp #{dirname}.tar.gz #{WEBSERVER}:#{FTP_DIR}/#{RELEASE_DIR}"

  if File.exists?("#{dirname}/CHANGELOG")
    sh "scp #{dirname}/CHANGELOG #{WEBSERVER}:#{FTP_DIR}"
  end
end

