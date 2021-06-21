from miscSupports import validate_path, open_setter, decode_line, flatten
from csvObject import write_csv


class VCFObject:
    def __init__(self, path):

        self._path = validate_path(path)
        self._zipped = (self._path.suffix == ".gz")

        # Extract the raw headers, and find the start byte of the file
        self.vcf_headers = []
        self.data_headers, self.start_byte = self._extract_headers()
        self._data_dict = {header: index for index, header in enumerate(self.data_headers)}

        # Determine how many format columns their are
        if "FORMAT" in list(self._data_dict.keys()):
            self._format_length = (len(self._data_dict) - 1) - self._data_dict["FORMAT"]
        else:
            self._format_length = 0

        # Extract known header titles as dicts
        self._info_dict = self._set_vcf_header("INFO=", "info")
        self._format_dict = self._set_vcf_header("FORMAT=", "format")
        self._meta_dict = self._set_vcf_header("META=", "meta")

        all_headers = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"] + list(self._info_dict.keys()) + \
                        flatten([[f"{header}_{index}" for header in list(self._format_dict.keys())]
                            for index in range(self._format_length)])

        self.header_dict = {header: index for index, header in enumerate(all_headers)}

    def _extract_headers(self):
        """
        Extract VCF and column headers from the VCF file

        :return: The headers of the columns and the start byte of the data
        """

        with open_setter(self._path)(self._path) as file:

            for line_byte in file:

                # Decode line
                line = decode_line(line_byte, self._zipped, "\t")

                # Then we have found the column headers and can stop
                if len(line) != 1:
                    return [self._clean(header) for header in line], file.tell()

                # Append the VCF header to the list of VCF headers
                else:
                    self.vcf_headers.append(self._clean(line[0]))

    @staticmethod
    def _clean(text):
        """
        Clean a string from vcf of characters that may be present

        Note
        ----
        VCF files uses ## and # to denote VCF headers and data headers respectively. Both \n and \t may be present
        """
        return text.replace("#", "").replace("\n", "").replace("\t", "")

    def _set_vcf_header(self, header_separator, header_unique):
        """
        For a given header separator extract a dict of {ID: {Number, Type, Description}}

        :param header_separator: String to separate on
        :type header_separator: str

        :param header_unique: Some headers may duplicate so this adds a value to it, this way each header will be unique
        :type header_unique: str

        :return: A dict of {ID: {Number, Type, Description}} for this separator
        :rtype: dict
        """

        header_dict = {}
        index = 0
        for header in self.vcf_headers:
            if header_separator in header:
                # Extract the line with the header value, then split the values on = to get the known output
                var_id, num, var_type, description = \
                    [v.split("=")[1].replace(">", "") for v in header.split(header_separator)[1].split(",")]

                header_dict[f"{var_id}_{header_unique}"] = \
                    {"Number": num, "Type": self._set_type(var_type), "Description": description, "Index": index}
                index += 1

        return header_dict

    @staticmethod
    def _set_type(variable_type):
        """Types may not be pass eval, so check via eval and then revert to backup comparision before failing"""

        try:
            return eval(variable_type.lower())
        except NameError:
            if variable_type.lower() == "string":
                return str
            elif variable_type.lower() == "integer":
                return int
            else:
                raise TypeError(f"Unexpected type of {variable_type} found in header")

    def covert_to_summary(self, info=True):

        out_rows = []
        with open_setter(self._path)(self._path) as file:

            file.seek(self.start_byte)
            for count_index, line_byte in enumerate(file):

                if count_index % 100000 == 0:
                    print(f"Extracted {count_index} lines")

                # Decode line
                line = decode_line(line_byte, self._zipped, "\t")

                row = line[:7] + ["NA" for _ in range(len(self._info_dict.keys()))] + \
                      ["NA" for _ in range(len(self._format_dict.keys()) * self._format_length)]

                if ("INFO" in self.data_headers) and info:

                    # Extract the info parameters
                    parameters = line[self._data_dict["INFO"]].split(";")

                    for para in parameters:
                        key, value = para.split("=")
                        row[self.header_dict[f"{key}_info"]] = value

                if "FORMAT" in self.data_headers:
                    parameter_names = line[self._data_dict["FORMAT"]].split(":")

                    for index, i in enumerate(range(self._data_dict["FORMAT"] + 1, len(line))):
                        parameter_values = line[i].replace("\n", "").split(":")

                        for name, value in zip(parameter_names, parameter_values):
                            row[self.header_dict[f"{name}_format_{index}"]] = value

                out_rows.append(row)

        write_csv(r"C:\Users\Samuel\Documents\Genetic_Examples", "Test", list(self.header_dict.keys()), out_rows)



if __name__ == '__main__':
    VCFObject(r"C:\Users\Samuel\Documents\Genetic_Examples\ukb-b-6134.vcf.gz").covert_to_summary()
