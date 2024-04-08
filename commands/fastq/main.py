
import argparse

import modules as m


class Main():

  @staticmethod
  def make_mock_file(fpath):

    fastq_file = m.Fastq(fpath)

    # file_size = m.file_sys.File.get_size(fastq_file.fpath, 'B')
    # print(f'file_size: {file_size}')
    # file_size = m.file_sys.File.get_size(fastq_file.fpath, 'KB')
    # print(f'file_size: {file_size}')
    # file_size = m.file_sys.File.get_size(fastq_file.fpath, 'MB')
    # print(f'file_size: {file_size}')
    file_size = m.file_sys.File.get_size(fastq_file.fpath, 'GB')
    print(f'file_size: {file_size} GB')
    # file_size = m.file_sys.File.get_size(fastq_file.fpath, 'TB')
    # print(f'file_size: {file_size}')

    # base_qty, seq_qty = fastq_file.get_base_and_seq_qty()
    # print(f'base_qty: {base_qty}')
    # print(f'base_qty: {base_qty/1_000_000_000} Gb')
    # print(f'seq_qty: {seq_qty}')

    mock_fastq_file = fastq_file.generate_mock_from_self(
      "mock-data-30Gbp.fastq",
      base_qty = 30 * 1_000_000_000
    )

    print(f'mock_fastq_file.fpath: {mock_fastq_file.fpath}')


if __name__ == "__main__":

  arg_parser = argparse.ArgumentParser()
  subparsers = arg_parser.add_subparsers()

  arg_parser.add_argument(
    'command_name',
    choices = [
      'mock'
    ]
  )
  args = arg_parser.parse_args()
  command_name = args.command_name

  # subparser_mock = subparsers.add_parser()
  # subparser_mock.add_argument(
  #   "-f",
  #   "--filepath",
  #   type = str
  # )
  # args_mock = subparser_mock.parse_args()
  

  match command_name:

    case 'mock':
      
      print(f'command_name: {command_name}')
      # fpath = args_mock.filepath

      # Main.make_mock_file(fpath)