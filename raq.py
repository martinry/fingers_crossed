from flask import Flask, redirect, render_template, request, session
from flask import url_for, send_from_directory
from collections import defaultdict
from pathlib import Path
import requests, sys
import json
import pyensembl
import subprocess
import os
import itertools

app = Flask(__name__)
class pfam_register():
    def __init__(self, 
                pfam_id,
                pfam_name,
                pfam_type):
        self.pfam_id = pfam_id
        self.pfam_name = pfam_name
        self.pfam_type = pfam_type
        self.positions = []
    
    def add_position(self, seq_from, seq_to, significance):
        new_position = [seq_from, seq_to, significance]
        if new_position not in self.positions:
            self.positions.append(new_position)
            
    def get_positions(self):
        return self.positions
    
class protein_register():
    def __init__(self,
                prot_id):
        self.prot_id = prot_id
        self.pfam_register_list = []
        
    def add_pfam_register(self, pfam_register):
        if pfam_register not in self.pfam_register_list:
            self.pfam_register_list.append(pfam_register)
        
    def get_id(self):
        return self.prot_id
    
    def get_pfam_registers_list(self):
        return self.pfam_register_list

def get_orthologs(source_id):
    server = "http://rest.ensembl.org"
    ext = "/homology/id/" + source_id + "?sequence=cdna;type=orthologues"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
      r.raise_for_status()
      sys.exit()
    return (r.json())

def get_proteins(seq_id):
    server = "http://rest.ensembl.org"
    ext = "/sequence/id/" + seq_id + "?multiple_sequences=1;type=protein"

    r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()
    encoded_proteins = []
    seq = ''
    prot_id = ''
    id_found = False
    
    proteins = r.text.split('>')[1:]
    ids = [i.split('\n')[0] for i in proteins]
    proteins = [''.join(i.split('\n')[1:]) for i in proteins]
    
    for i in range(0, len(ids)):
        encoded_proteins.append([ids[i], proteins[i]])
    
    file_name = seq_id + '.fasta'
    fout = open(file_name, 'w')
    fout.write(r.text)
    
    return encoded_proteins

# make app 

def parse_pfam(gene_id):
    prot_objects = {}
    
    fin = gene_id + '.pfam'
    with open(fin) as pfam:
        for line in pfam:
            if line[0] != '#' and line != '\n':
                line = '\t'.join(line.split())
                line_array = line.split('\t')
                
                prot_id = line_array[0]
                seq_from = line_array[1]
                seq_to = line_array[2]
                pfam_id = line_array[5]
                pfam_name = line_array[6]
                pfam_type = line_array[7]
                significance = line_array[12]
                
                if prot_id not in list(prot_objects.keys()):
                    prot_objects[prot_id] = protein_register(prot_id)
                
                current_pfams_ids = [i.pfam_id for i in prot_objects[prot_id].get_pfam_registers_list()]
                
                if pfam_id not in current_pfams_ids:
                    pfam_object = pfam_register(pfam_id, pfam_name, pfam_type)
                    pfam_object.add_position(seq_from, seq_to, significance) 
                    prot_objects[prot_id].add_pfam_register(pfam_object)
                else:
                    for i in prot_objects[prot_id].get_pfam_registers_list():
                        if i.pfam_id == pfam_id:
                            i.add_position(seq_from, seq_to, significance)
                            break
                    
    return(prot_objects)    

alignment = ''
@app.route('/')
def index():
    return render_template('index.html')


@app.route('/compare', methods=['GET', 'POST'])
def compare():
    if request.method == 'GET':
        return render_template('compare.html')
    else:
        form = request.form
        source_id = form['source_id']
        ortho_id = form['ortho_id']
        
        if (pyensembl.common.is_valid_ensembl_id(source_id)) and (pyensembl.common.is_valid_ensembl_id(ortho_id)): # good user
            
            json_data = get_orthologs(source_id)
            
           # o_json = open('ortho.json', 'w')
            #o_json.write(str(json_data))
            #o_json.close()
            
            try:
                orthologs = [i['target']['id'] for i in json_data['data'][0]['homologies']]    
            except:
                error = "I'm sorry, ensembl does not provide a list of orthologs for the id " + source_id + '.'
                return render_template('input_error.html', error = error)
            
            if ortho_id in orthologs:
                info = [i for i in json_data['data'][0]['homologies'] if i['target']['id'] == ortho_id] 
                ortho = info[0]['target']
                ortho_gene = ortho['id']
                ortho_taxon = ortho['taxon_id']
                ortho_cigar = ortho['cigar_line']
                ortho_specie = ortho['species'].replace('_', ' ').title()
                ortho_alignment = ortho['align_seq']
                
                source = info[0]['source']
                source_gene = source_id
                source_taxon = source['taxon_id']
                source_specie = source['species'].replace('_', ' ').title()
                source_alignment = source['align_seq']
                
                perc_id = ortho['perc_id']
                
                a = 0
                total_alignment = []
                for i in range (0, max(len(ortho_alignment), len(source_alignment))):
                    if i%90 == 0 or i == max(len(ortho_alignment), len(source_alignment)):
                        total_alignment.append([ortho_alignment[a:i], source_alignment[a:i]])
                        a = i
                total_alignment = total_alignment[1:]
                
                #print('\n\n\n\n\n\n\n' + os.getcwd() + '\n\n\n\n\n\n')
                directory = os.path.join(os.getcwd(), 'data')
                os.chdir(directory)
                
                ortho_proteins = get_proteins(ortho_id)
                source_proteins = get_proteins(source_id)     
                
                import paramiko
                from scp import SCPClient
                
                files = [f for f in os.listdir(directory) if f.endswith('fasta')]
                
                pfam_files = [f.replace('fasta', 'pfam') for f in files]
                for f in pfam_files:
                    print(f)
                    if Path(f).is_file():
                        pass
                    else:
                        ssh = paramiko.SSHClient()
                        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
                        ssh.connect('130.235.244.18', username='ms1705', password='GaTaCA1705')
                        scp = SCPClient(ssh.get_transport())

                        print(f)
                        id_file = f.replace('.pfam', '')
                        print(id_file)
                        fasta_f = f.replace('.pfam', '.fasta')
                        scp.put(fasta_f, 'pfam_test')
                        print('On server!')

                        command = "pfam_scan.pl -fasta ./pfam_test/" + fasta_f + " -outfile ./pfam_test/" + id_file + ".pfam -e_seq 1e-5 -e_dom 1e-5 -dir /opt/bio/PFAMdb"
                        print(command)
                        ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command(command)
                        print(ssh_stderr.read().decode('utf-8'))
                        while ssh_stdout.channel.recv_exit_status() != 0:
                            None

                        sftp = ssh.open_sftp()
                        file_remote = './pfam_test/' + id_file + '.pfam'
                        file_local = id_file + '.pfam'
                        sftp.get(file_remote, file_local)
                        
                        files_in_remote = sftp.listdir(path='./pfam_test/')
                        for file in files_in_remote:
                            sftp.remove('./pfam_test/'+file)
                        
                        print('Done!')

                        sftp.close()
                        ssh.close()
                
                ortho_pfam = parse_pfam(ortho_id)
                source_pfam = parse_pfam(source_id)
                os.chdir('..')
                
                return render_template('ortho_info.html', total_align = total_alignment, source_gene = source_gene, ortho_gene = ortho_gene, source_taxon = source_taxon, ortho_taxon = ortho_taxon, ortho_cigar = ortho_cigar, source_specie = source_specie, ortho_specie = ortho_specie, perc_id = perc_id, ortho_proteins = ortho_proteins, source_proteins = source_proteins, ortho_pfam = ortho_pfam, source_pfam = source_pfam)
            else:
                error = 'Sorry, these ids are not registered as orthologs at ensembl.' + "\n" + 'We found these orthologs for the id ' + source_id + ': ' + str(orthologs) 
                return render_template('input_error.html', error = error)
        
        elif (pyensembl.common.is_valid_ensembl_id(source_id) == False) and (pyensembl.common.is_valid_ensembl_id(ortho_id) == False):
            error = "I'm sorry, the ids you entered doesn't seem to be valid ensembl ids."
            return render_template('input_error.html', error = error) # both wrong
        
        elif (pyensembl.common.is_valid_ensembl_id(ortho_id)):
            error = "I'm sorry, the id " + source_id + " doesn't seem to be a valid ensembl id."
            return render_template('input_error.html', error = error) # your first id is wrong
        
        else:
            error = "I'm sorry, the id " + ortho_id + " doesn't seem to be a valid ensembl id."
            return render_template('input_error.html', error = error) # your second id is wrong


if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5003, debug=True)