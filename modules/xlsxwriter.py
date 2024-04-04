
import xlsxwriter 

import config
import modules as m


class XlsxWriter():
    

  @staticmethod
  def tsv_to_xlsx(tsv_filename, xlsx_filename):
    # Create a new Excel workbook and add a worksheet
    workbook = xlsxwriter.Workbook(xlsx_filename)
    worksheet = workbook.add_worksheet("Bacfile")

    # Open the TSV file for reading
    with open(tsv_filename, 'r') as tsv_file:
      # Read each line from the TSV file
      for row_num, line in enumerate(tsv_file):
        # Split the line into fields using the tab character
        fields = line.strip().split('\t')

        # Write each field to the worksheet
        for col_num, field in enumerate(fields):

          field = m.String.try_str_to_number(field)

          worksheet.write(row_num, col_num, field)

    # Close the workbook
    workbook.close()