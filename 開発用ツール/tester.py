# テスト用ツール
# ./clay_tester の中にある hoge.cpp を実行してエラーなく終了するか試す
# hoge.in （というファイルがあれば）を標準入力から食わせ，標準出力が hoge.out と一致しないとエラーとする
# (hoge.out が存在しなければ，標準出力は空になるはずとみなす)

import datetime
import os
import random, string
import subprocess
import glob
from natsort import natsorted
from concurrent.futures import ProcessPoolExecutor

def run_command(command):
  try:
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True, timeout=60)
  except subprocess.TimeoutExpired as e:
    result = subprocess.run('echo a', stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True, timeout=60)
    result.returncode = 1
    result.stdout = 'timeout'
    result.stderr = 'timeout'
  return result

def randomname(n):
  randlst = [random.choice(string.ascii_letters + string.digits) for i in range(n)]
  return ''.join(randlst)

def file_add_writer(filename, wstr):
  dt_now = datetime.datetime.now()
  dt_str = dt_now.strftime('%Y%m%d-%H%M%S')
  with open(filename, 'a') as f:
    print(dt_str + ': ' + wstr, file=f)

def file_reader(filename):
  f = open(filename, 'r')
  res = f.read()
  f.close()
  return res

def test_doit(filename):
  print('test_begin: ' + filename)
  exefilename = 'test_' + randomname(10)
  
  com = './clayexe.exe < clay_tester/' + filename + '.cpp > ' + exefilename + '.cpp'
  res = run_command(com)

  file_add_writer('tester.log', 'run: '+com)
  file_add_writer('tester.log', 'returncode: '+ str(res.returncode))

  if res.returncode != 0:
    file_add_writer('tester_failed.log', 'run: ' + com)
    file_add_writer('tester_failed.log', 'returncode: '+ str(res.returncode))
    file_add_writer('tester_failed.log', 'stdout: ' + res.stdout)
    file_add_writer('tester_failed.log', 'stderr: ' + res.stderr)
    file_add_writer('tester_failed_list.log', filename)
    return
  
  com = 'g++ -O2 -std=gnu++17 -o ' + exefilename + '.exe ' + exefilename + '.cpp'
  res = run_command(com)

  file_add_writer('tester.log', 'run: '+com)
  file_add_writer('tester.log', 'returncode: '+ str(res.returncode))

  com = 'rm ' + exefilename + '.cpp'
  run_command(com)

  if res.returncode != 0:
    file_add_writer('tester_failed.log', 'run: ' + com)
    file_add_writer('tester_failed.log', 'returncode: '+ str(res.returncode))
    file_add_writer('tester_failed.log', 'stdout: ' + res.stdout)
    file_add_writer('tester_failed.log', 'stderr: ' + res.stderr)
    file_add_writer('tester_failed_list.log', filename)
    return

  com = './' + exefilename + '.exe'
  if os.path.isfile('./clay_tester/'+filename+'.in'):
    com = com + ' < ' + './clay_tester/' + filename + '.in'

  res = run_command(com)
  file_add_writer('tester.log', 'run: '+com)
  file_add_writer('tester.log', 'returncode: '+ str(res.returncode))

  com = 'rm ' + exefilename + '.exe'
  run_command(com)

  out_test = res.stdout
  out_corr = ''
  if os.path.isfile('./clay_tester/'+filename+'.out'):
    out_corr = file_reader('./clay_tester/'+filename+'.out')

  if res.returncode != 0 or out_test != out_corr:
    file_add_writer('tester_failed.log', 'run: ' + com)
    file_add_writer('tester_failed.log', 'returncode: '+ str(res.returncode))
    file_add_writer('tester_failed.log', 'stdout: ' + res.stdout)
    file_add_writer('tester_failed.log', 'stderr: ' + res.stderr)
    file_add_writer('tester_failed.log', 'correct_out: ' + out_corr)
    file_add_writer('tester_failed_list.log', filename)
    return


  file_add_writer('tester_success_list.log', filename)



files = glob.glob("./clay_tester/*.cpp")
for i in range(len(files)):
  files[i] = files[i].replace('./clay_tester/', '')
  files[i] = files[i].replace('.cpp', '')

files = natsorted(files)

with ProcessPoolExecutor(max_workers=4) as executor:
  results = list(executor.map(test_doit, files))

