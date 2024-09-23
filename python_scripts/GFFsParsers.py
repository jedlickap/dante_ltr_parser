class GFF3Line:
    def __init__(self, line):
        self.raw_line = line
        if ";" == self.raw_line.rstrip()[-1]: 
            self.raw_line = self.raw_line[:-2]
        self.fields = self.raw_line.split('\t')
        self.seqid = self.fields[0]
        self.source = self.fields[1]
        self.type = self.fields[2]
        self.start = int(self.fields[3])
        self.end = int(self.fields[4])
        self.score = self.fields[5]
        self.strand = self.fields[6]
        self.phase = self.fields[7]
        self.attributes = self._parse_attributes()

    def _parse_attributes(self):
        attributes_str = self.fields[8]
        attributes = {}

        if attributes_str != '.':
            key_value_pairs = attributes_str.split(';')
            for pair in [pr for pr in key_value_pairs if pr]:
                key, value = pair.split('=')
                attributes[key] = value

        return attributes

    def get_attribute(self, key):
        return self.attributes.get(key)

class DanteLTrGFF:
    
    def __init__(self, gff3_path):
        self.file_path = gff3_path
        self.gff_lines = self._parse_gff_file()

    def _parse_gff_file(self):
        parsed_lines = []
        with open(self.file_path, 'r') as file:
            for line in file:
                if not line.startswith("#"):
                    gff_line = GFF3Line(line)
                    parsed_lines.append(gff_line)
        return parsed_lines

    def _get_id(self, gff_line):
        if gff_line.get_attribute('ID'): 
                return gff_line.get_attribute('ID').strip()
        else:
                return gff_line.get_attribute('Parent').strip()
        
    def process_data(self):
        tes_dict = {}
        for gff_line in self.gff_lines:
            # Do something with each GFF line, e.g., print attributes
            TE_ID = self._get_id(gff_line)
            if TE_ID not in tes_dict: 
                tes_dict[TE_ID] = {'transposable_element':[],
                                'long_terminal_repeat':[],
                                'protein_domain':[],
                                'primer_binding_site':[],
                                'target_site_duplication':[]}
                tes_dict[TE_ID][gff_line.type].append(gff_line)
            else:
                tes_dict[TE_ID][gff_line.type].append(gff_line)
        return tes_dict
        
# Example usage
    # dante_instance = DanteLTrGFF('DANTE_LTR_out.gff3')
    # tes_dict = dante_instance.process_data()
    # print(tes_dict['TE_00000001_CgrY'])
    # for ID in tes_dict:
    #     if tes_dict[ID]['transposable_element']:
    #         gff_line = tes_dict[ID]['transposable_element'][0]
    #         if 'T' in gff_line.attributes['Rank']:
    #             print(f"{ID}\t{gff_line.attributes['Rank']}")
            
# print(tes_dict)

class EdtaGFF:
    
    def __init__(self, gff3_path):
        self.file_path = gff3_path
        self.gff_lines = self._parse_gff_file()

    def _parse_gff_file(self):
        parsed_lines = []
        with open(self.file_path, 'r') as file:
            for line in file:
                if not line.startswith("#"):
                    gff_line = GFF3Line(line)
                    parsed_lines.append(gff_line)
        return parsed_lines
        
    def process_data(self):
        """
        Keys in EDTA attriburtes:
        ID,Name,Classification,Sequence_ontology,Identity,Method
        ID,Name,Classification,Sequence_ontology,ltr_identity,Method
        ID,Parent,Name,Classification,Sequence_ontology,ltr_identity
        """
        tes_dict = {}
        for gff_line in self.gff_lines:
            if gff_line.type not in tes_dict: 
                tes_dict[gff_line.type] = []
                tes_dict[gff_line.type].append(gff_line)
            else:
                tes_dict[gff_line.type].append(gff_line)
        return tes_dict
