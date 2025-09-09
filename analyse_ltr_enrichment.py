#!/usr/bin/env python3
import sys
import os
from collections import defaultdict
import re
import random 
import shutil 
import subprocess

# Sprawdzenie i obsługa importów na początku
try:
    import pyfastx
except ImportError:
    print("BŁĄD: Biblioteka pyfastx nie jest zainstalowana. Jest wymagana do wczytywania FASTA.")
    print("Zainstaluj ją używając: pip install pyfastx")
    sys.exit(1)

try:
    from scipy.stats import fisher_exact
except ImportError:
    print("BŁĄD: Biblioteka SciPy (scipy) nie jest zainstalowana. Jest wymagana do testu Fishera.")
    print("Zainstaluj ją używając: pip install scipy")
    sys.exit(1)

try:
    from translate_seq import six_frame_translate 
except ImportError:
    # Ten import jest potrzebny tylko jeśli Etap 1 jest wykonywany.
    # Jeśli Etap 1 jest pomijany, brak tego modułu nie powinien być krytyczny,
    # ale dla kompletności skryptu zostawiamy go z ostrzeżeniem.
    print("OSTRZEŻENIE: Nie można zaimportować six_frame_translate z translate_seq.py.")
    print("Jeśli Etap 1 (analiza HMMER) ma być wykonany, ten plik jest niezbędny.")
    # Nie robimy sys.exit(1) tutaj, aby umożliwić pominięcie Etapu 1.
    # Funkcja six_frame_translate będzie wołana tylko w run_hmmsearch_on_ltrs.

SCRIPT_PATH_BASE = os.path.dirname(os.path.abspath(__file__))

def parse_gff_attributes(attributes_str):
    attributes = {}
    for part in attributes_str.split(';'):
        if '=' in part:
            key, value = part.split('=', 1)
            attributes[key.strip()] = value.strip()
    return attributes

def get_chromosome_lengths(fasta_file):
    lengths = {}
    try:
        fa = pyfastx.Fasta(fasta_file)
        for record in fa:
            lengths[record.name] = len(record.seq)
    except Exception as e:
        print(f"BŁĄD podczas wczytywania pliku FASTA {fasta_file}: {e}")
        sys.exit(1) 
    
    if not lengths: 
        print(f"OSTRZEŻENIE: Nie wczytano żadnych sekwencji (lub długości) z pliku FASTA: {fasta_file}")
        print("Upewnij się, że plik FASTA jest poprawny i nie jest pusty.")
        sys.exit(1)
    return lengths

def calculate_overlap(start1, end1, start2, end2):
    return max(0, min(end1, end2) - max(start1, start2) + 1)

def run_hmmsearch_on_ltrs(input_ltr_fasta_file, output_hmm_tblout_file, script_base_path):
    print(f"Rozpoczynanie translacji i analizy HMMER dla pliku: {input_ltr_fasta_file}")
    
    translated_ltr_file = input_ltr_fasta_file + ".translated.fasta"
    
    if not os.path.exists(input_ltr_fasta_file) or os.path.getsize(input_ltr_fasta_file) == 0:
        print(f"OSTRZEŻENIE: Plik z sekwencjami LTR ({input_ltr_fasta_file}) jest pusty lub nie istnieje. Pomijanie translacji i HMMER.")
        open(translated_ltr_file, 'w').close() 
        open(output_hmm_tblout_file, 'w').close()
        return False 

    # Sprawdzenie czy translate_seq jest dostępne, zanim go użyjemy
    if 'six_frame_translate' not in globals():
        print("BŁĄD KRYTYCZNY: Funkcja six_frame_translate nie jest dostępna (problem z importem translate_seq.py).")
        return False

    try:
        with open(translated_ltr_file, "w") as f_trans_out:
            six_frame_translate(input_ltr_fasta_file, f_trans_out)
        print(f"Zakończono translację do: {translated_ltr_file}")
    except Exception as e_translate:
        print(f"Błąd podczas translacji sekwencji LTR z {input_ltr_fasta_file}: {e_translate}")
        return False

    if not os.path.exists(translated_ltr_file) or os.path.getsize(translated_ltr_file) == 0:
        print(f"OSTRZEŻENIE: Plik z przetłumaczonymi sekwencjami LTR ({translated_ltr_file}) jest pusty. Pomijanie hmmsearch.")
        open(output_hmm_tblout_file, 'w').close()
        return False

    hmmsearch_executable = os.path.join(script_base_path, 'bin', 'hmmsearch')
    rexdb_hmm_file = os.path.join(script_base_path, 'bin', 'REXdb.hmm')

    if not os.path.exists(hmmsearch_executable):
        print(f"BŁĄD: Nie znaleziono programu hmmsearch w {hmmsearch_executable}")
        return False
    if not os.path.exists(rexdb_hmm_file):
        print(f"BŁĄD: Nie znaleziono pliku bazy danych REXdb.hmm w {rexdb_hmm_file}")
        return False

    hmm_cmd = [
        hmmsearch_executable, '--noali', '--tblout', output_hmm_tblout_file,
        '-E', '1e-5', '--cpu', '4', 
        rexdb_hmm_file, translated_ltr_file
    ]
    print(f"Uruchamianie HMMER: {' '.join(hmm_cmd)}")
    try:
        proc_hmm = subprocess.Popen(hmm_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout_hmm, stderr_hmm = proc_hmm.communicate()
        if proc_hmm.returncode != 0:
            print(f"Błąd podczas uruchamiania HMMER:\n{stderr_hmm.decode(errors='replace')}")
    except Exception as e_hmm:
        print(f"Wyjątek podczas uruchamiania HMMER: {e_hmm}")
        return False
    
    print(f"Zakończono analizę HMMER. Surowe wyniki w: {output_hmm_tblout_file}")
    return True


def process_hmmsearch_results(hmm_tblout_file, output_ltr_domain_annotations_file):
    print(f"Przetwarzanie wyników HMMER z {hmm_tblout_file} do {output_ltr_domain_annotations_file}")
    processed_hits = set() 
    
    if not os.path.exists(hmm_tblout_file):
        print(f"OSTRZEŻENIE: Plik wyników HMMER {hmm_tblout_file} nie został utworzony lub jest pusty. Plik adnotacji domen będzie pusty.")
        open(output_ltr_domain_annotations_file, 'w').close()
        return

    with open(hmm_tblout_file) as f_hmm_in, open(output_ltr_domain_annotations_file, "w") as f_final_out:
        f_final_out.write("LTR_ID\tDomain_Type\n") 
        for line in f_hmm_in:
            if not line.startswith("#"):
                parts = line.strip().split()
                if len(parts) >= 3: 
                    query_name_hmm = parts[0] 
                    domain_name_hmm = parts[2] 
                    ltr_id_from_fasta_header = query_name_hmm.split("|")[0]
                    hit_key = f"{ltr_id_from_fasta_header}\t{domain_name_hmm}"
                    if hit_key not in processed_hits:
                         f_final_out.write(hit_key + "\n")
                         processed_hits.add(hit_key)
    print(f"Zapisano przetworzone adnotacje domen do: {output_ltr_domain_annotations_file}")


def main():
    try:
        current_script_path = os.path.abspath(__file__)
        base_dir = os.path.dirname(current_script_path)
        print(f"INFO: Automatycznie wykryty base_dir (katalog skryptu): {base_dir}")
        if not os.path.isdir(os.path.join(base_dir, "bin")):
             print(f"OSTRZEŻENIE: W automatycznie wykrytym base_dir ({base_dir}) nie ma podkatalogu 'bin'.")
             print("Jeśli skrypt nie znajduje się w głównym katalogu CentIER, zmodyfikuj 'base_dir' ręcznie.")
             # Próba użycia katalogu roboczego jako alternatywy, jeśli 'bin' tam jest
             if os.path.isdir(os.path.join(os.getcwd(), "bin")):
                 base_dir = os.getcwd()
                 print(f"INFO: Zmieniono base_dir na bieżący katalog roboczy: {base_dir}")
             else: # Ostateczny fallback
                 base_dir = "/media/kaloryfer/6cd6dcf9-ce32-4bda-8b85-1f45b4f85630/Pellialong/annotation/egapx_output/CentIER/"
                 print(f"OSTRZEŻENIE: Używam domyślnej (hardkodowanej) ścieżki base_dir: {base_dir}")
    except NameError: 
        base_dir = "/media/kaloryfer/6cd6dcf9-ce32-4bda-8b85-1f45b4f85630/Pellialong/annotation/egapx_output/CentIER/"
        print(f"OSTRZEŻENIE: Nie udało się automatycznie wykryć base_dir. Używam domyślnej: {base_dir}")
    print("INFO: Upewnij się, że ścieżka base_dir jest poprawnym katalogiem głównym CentIER.")
    print("INFO: Jeśli nie, zmodyfikuj zmienną 'base_dir' w skrypcie lub uruchom skrypt z katalogu CentIER.\n")
    
    genome_fasta_file = os.path.join(base_dir, "complete.genomic.fna")
    centromere_ranges_file = os.path.join(base_dir, "CentIER_final_results", "complete.genomic.fna_centromere_range.txt")
    genome_ltr_gff3_file = os.path.join(base_dir, "complete.genomic.fna.LTR.gff3") 
    
    all_genome_ltrs_fasta_file = os.path.join(base_dir, "CentIER_final_results", "complete.genomic.fna.all_genome_ltrs.fasta")
    all_genome_ltrs_hmm_tblout_file = os.path.join(base_dir, "CentIER_final_results", "complete.genomic.fna.all_genome_ltrs.hmm.tblout")
    all_genome_ltr_domain_annotations_file = os.path.join(base_dir, "CentIER_final_results", "complete.genomic.fna.all_genome_ltr_rexdb_domains.tsv")
    output_report_file = os.path.join(base_dir, "CentIER_final_results", "genome_wide_rexdb_domain_enrichment_report.tsv") 
    
    min_ltr_len_param = 50
    p_value_threshold_param = 0.05
    top_n_to_report = 10 
    
    # --- KLUCZOWY PARAMETR DO KONTROLOWANIA ETAPU 1 ---
    # Ustaw na False, jeśli chcesz wymusić ponowne wykonanie Etapu 1 nawet jeśli pliki istnieją
    skip_stage1_if_files_exist = True 
    # --- ---

    for f_path_check in [genome_fasta_file, centromere_ranges_file, genome_ltr_gff3_file]:
        if not os.path.exists(f_path_check):
            print(f"BŁĄD KRYTYCZNY: Niezbędny plik wejściowy nie istnieje: {f_path_check}")
            sys.exit(1)
            
    print(f"Wczytywanie długości chromosomów z: {genome_fasta_file}")
    chr_lengths = get_chromosome_lengths(genome_fasta_file)
    total_genome_length = sum(chr_lengths.values())
    print(f"Całkowita długość genomu: {total_genome_length} bp")

    print(f"Wczytywanie regionów centromerowych z: {centromere_ranges_file}")
    centromere_regions = defaultdict(list) 
    total_centromeric_length = 0
    try:
        with open(centromere_ranges_file, 'r') as f_cen:
            for line in f_cen:
                if line.strip() and not line.startswith("#") and not line.startswith("[source"):
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        chrom, start_str, end_str = parts[0], parts[1], parts[2]
                        try:
                            start, end = int(start_str), int(end_str)
                            if start > end : start, end = end, start 
                            centromere_regions[chrom].append((start, end)) 
                            total_centromeric_length += (end - start + 1)
                        except ValueError:
                            print(f"Ostrzeżenie: Pomijanie linii z niepoprawnymi koordynatami w pliku regionów centromerowych: {line.strip()}")
                    else:
                        print(f"Ostrzeżenie: Pomijanie niepoprawnej linii w pliku regionów centromerowych: {line.strip()}")
    except FileNotFoundError:
        print(f"BŁĄD: Nie znaleziono pliku regionów centromerowych: {centromere_ranges_file}")
        sys.exit(1)
    
    if total_centromeric_length == 0:
        print(f"BŁĄD: Brak zdefiniowanych regionów centromerowych lub plik {centromere_ranges_file} jest pusty/niepoprawny.")
        sys.exit(1)
        
    total_non_centromeric_length = total_genome_length - total_centromeric_length
    print(f"Całkowita długość regionów centromerowych: {total_centromeric_length} bp ({total_centromeric_length/total_genome_length*100:.2f}% genomu)")
    print(f"Całkowita długość regionów niecentromerowych: {total_non_centromeric_length} bp ({total_non_centromeric_length/total_genome_length*100:.2f}% genomu)")

    print("\n--- ETAP 1: Generowanie całogenomowej adnotacji domen REXdb dla LTR ---")
    if skip_stage1_if_files_exist and os.path.exists(all_genome_ltr_domain_annotations_file) and os.path.getsize(all_genome_ltr_domain_annotations_file) > len("LTR_ID\tDomain_Type\n"): # Sprawdź czy nie tylko nagłówek
        print(f"INFO: Plik {all_genome_ltr_domain_annotations_file} już istnieje i zawiera dane. Pomijanie Etapu 1.")
    else:
        if 'six_frame_translate' not in globals():
            print("BŁĄD KRYTYCZNY: Funkcja six_frame_translate nie jest zaimportowana (brak translate_seq.py?). Nie można wykonać Etapu 1.")
            sys.exit(1)

        print(f"INFO: Rozpoczynanie Etapu 1: Ekstrakcja LTR, translacja i analiza HMMER.")
        genome_fa = pyfastx.Fasta(genome_fasta_file) 

        print(f"Ekstrahowanie sekwencji LTR z {genome_ltr_gff3_file} do {all_genome_ltrs_fasta_file}...")
        ltr_count_for_fasta = 0
        with open(genome_ltr_gff3_file, 'r') as f_gff, open(all_genome_ltrs_fasta_file, 'w') as f_fasta_out:
            for line_num, line in enumerate(f_gff):
                if line.startswith("#") or not line.strip(): continue
                fields = line.strip().split('\t')
                if len(fields) < 9: continue

                chrom, _, feature_type, start_str, end_str, _, strand, _, attributes_str = fields
                if feature_type not in ["LTR_retrotransposon", "Gypsy_LTR_retrotransposon", "Copia_LTR_retrotransposon", "repeat_region"]:
                    continue
                
                try: start, end = int(start_str), int(end_str)
                except ValueError: continue
                
                if (end - start + 1) < min_ltr_len_param: continue

                attributes = parse_gff_attributes(attributes_str)
                ltr_id = attributes.get("ID", f"LTR_{chrom}_{start}_{end}") 
                ltr_id = re.sub(r'[^\w.-]', '_', ltr_id)

                try:
                    seq_data = genome_fa.fetch(chrom, (start, end)) 
                    if strand == '-':
                        seq_data = pyfastx.reverse_complement(seq_data)
                    
                    f_fasta_out.write(f">{ltr_id}\n{seq_data}\n")
                    ltr_count_for_fasta +=1
                except KeyError:
                    print(f"Ostrzeżenie: Nie znaleziono chromosomu '{chrom}' w pliku FASTA dla LTR {ltr_id}.")
                except Exception as e_fetch_gff:
                    print(f"Błąd podczas pobierania sekwencji dla LTR {ltr_id} ({chrom}:{start}-{end}): {e_fetch_gff}")

        print(f"Zapisano {ltr_count_for_fasta} sekwencji LTR do {all_genome_ltrs_fasta_file}")

        if ltr_count_for_fasta > 0:
            hmm_success = run_hmmsearch_on_ltrs(all_genome_ltrs_fasta_file, all_genome_ltrs_hmm_tblout_file, SCRIPT_PATH_BASE)
            if hmm_success:
                process_hmmsearch_results(all_genome_ltrs_hmm_tblout_file, all_genome_ltr_domain_annotations_file)
            else:
                print("BŁĄD: Analiza HMMER nie powiodła się. Nie można kontynuować z analizą wzbogacenia domen.")
                sys.exit(1)
        else:
            print("OSTRZEŻENIE: Nie znaleziono żadnych LTR-ów spełniających kryteria w pliku GFF. Plik adnotacji domen będzie pusty.")
            open(all_genome_ltr_domain_annotations_file, 'w').write("LTR_ID\tDomain_Type\n")

    print("\n--- ETAP 2: Analiza wzbogacenia domen REXdb w centromerach ---")
    
    if not os.path.exists(all_genome_ltr_domain_annotations_file) or os.path.getsize(all_genome_ltr_domain_annotations_file) <= len("LTR_ID\tDomain_Type\n"):
        print(f"BŁĄD: Plik z całogenomowymi adnotacjami domen LTR ({all_genome_ltr_domain_annotations_file}) jest pusty lub nie istnieje.")
        print("Nie można przeprowadzić analizy wzbogacenia.")
        sys.exit(1)

    ltr_coordinates_and_length = {} 
    print(f"Ponowne wczytywanie GFF dla mapowania ID LTR na współrzędne: {genome_ltr_gff3_file}")
    try:
        with open(genome_ltr_gff3_file, 'r') as f_gff:
            for line in f_gff:
                if line.startswith("#") or not line.strip(): continue
                fields = line.strip().split('\t')
                if len(fields) < 9: continue
                chrom, _, feature_type, start_str, end_str, _, _, _, attributes_str = fields
                if feature_type not in ["LTR_retrotransposon", "Gypsy_LTR_retrotransposon", "Copia_LTR_retrotransposon", "repeat_region"]:
                    continue
                try: start, end = int(start_str), int(end_str)
                except ValueError: continue
                
                ltr_len_val = end - start + 1 
                if ltr_len_val < min_ltr_len_param: continue

                attributes = parse_gff_attributes(attributes_str)
                ltr_id_from_gff = attributes.get("ID", f"LTR_{chrom}_{start}_{end}")
                ltr_id_from_gff = re.sub(r'[^\w.-]', '_', ltr_id_from_gff) 
                ltr_coordinates_and_length[ltr_id_from_gff] = (chrom, start, end, ltr_len_val)
    except FileNotFoundError:
        print(f"Krytyczny błąd: Nie znaleziono pliku GFF {genome_ltr_gff3_file} przy drugiej próbie odczytu.")
        sys.exit(1)
    
    print(f"Wczytano współrzędne dla {len(ltr_coordinates_and_length)} LTR-ów z pliku GFF.")

    print(f"DEBUG: Typ zmiennej 'centromere_regions' tuż przed pętlą to: {type(centromere_regions)}")
    ch1_example_present = "Obecny" if 'Ch1' in centromere_regions and centromere_regions['Ch1'] else "Brak lub pusty dla Ch1"
    print(f"DEBUG: Czy 'centromere_regions' zawiera dane (przykład dla Ch1)? {ch1_example_present}. Liczba chromosomów w dict: {len(centromere_regions)}")

    domain_data = defaultdict(lambda: {
        "len_cen": 0, "len_noncen": 0,
        "count_cen": 0, "count_noncen": 0,
        "ltr_ids_cen": set(), "ltr_ids_noncen": set() 
    })
    all_rexdb_domain_types_found = set()
    
    total_ltr_elements_with_domains_in_centromeres = set()
    total_ltr_elements_with_domains_outside_centromeres = set()
    total_ltr_bp_with_domains_in_centromeres = 0
    total_ltr_bp_with_domains_outside_centromeres = 0

    print(f"Wczytywanie całogenomowych adnotacji domen LTR z: {all_genome_ltr_domain_annotations_file}")
    with open(all_genome_ltr_domain_annotations_file, 'r') as f_domains:
        next(f_domains) 
        for line in f_domains:
            if not line.strip(): continue
            parts = line.strip().split('\t')
            if len(parts) < 2: continue
            ltr_id, domain_type = parts[0], parts[1]
            
            all_rexdb_domain_types_found.add(domain_type)

            if ltr_id not in ltr_coordinates_and_length:
                continue
            
            chrom, ltr_start, ltr_end, ltr_len_val_gff = ltr_coordinates_and_length[ltr_id] 
            
            is_in_centromere_flag = False
            accumulated_length_in_centromere_for_ltr = 0
            if chrom in centromere_regions: 
                for cen_start, cen_end in centromere_regions[chrom]:
                    overlap = calculate_overlap(ltr_start, ltr_end, cen_start, cen_end)
                    if overlap > 0:
                        is_in_centromere_flag = True
                        accumulated_length_in_centromere_for_ltr += overlap
            
            if is_in_centromere_flag:
                if ltr_id not in domain_data[domain_type]["ltr_ids_cen"]: 
                    domain_data[domain_type]["len_cen"] += ltr_len_val_gff 
                    domain_data[domain_type]["count_cen"] += 1
                    total_ltr_elements_with_domains_in_centromeres.add(ltr_id) 
                domain_data[domain_type]["ltr_ids_cen"].add(ltr_id) 
            else:
                if ltr_id not in domain_data[domain_type]["ltr_ids_noncen"]:
                    domain_data[domain_type]["len_noncen"] += ltr_len_val_gff
                    domain_data[domain_type]["count_noncen"] += 1
                    total_ltr_elements_with_domains_outside_centromeres.add(ltr_id)
                domain_data[domain_type]["ltr_ids_noncen"].add(ltr_id)
    
    for ltr_id_calc in total_ltr_elements_with_domains_in_centromeres: # Użyj innej nazwy zmiennej
        _, _, _, ltr_len_val_calc = ltr_coordinates_and_length[ltr_id_calc]
        total_ltr_bp_with_domains_in_centromeres += ltr_len_val_calc
    for ltr_id_calc in total_ltr_elements_with_domains_outside_centromeres: # Użyj innej nazwy zmiennej
        _, _, _, ltr_len_val_calc = ltr_coordinates_and_length[ltr_id_calc]
        total_ltr_bp_with_domains_outside_centromeres += ltr_len_val_calc

    print(f"Znaleziono {len(all_rexdb_domain_types_found)} unikalnych typów domen REXdb w LTR-ach.")
    print(f"Sumaryczna długość LTR zawierających domeny REXdb w centromerach: {total_ltr_bp_with_domains_in_centromeres} bp")
    print(f"Sumaryczna długość LTR zawierających domeny REXdb poza centromerami: {total_ltr_bp_with_domains_outside_centromeres} bp")
    print(f"Sumaryczna liczba unikalnych LTR-ów z domenami REXdb w centromerach: {len(total_ltr_elements_with_domains_in_centromeres)}")
    print(f"Sumaryczna liczba unikalnych LTR-ów z domenami REXdb poza centromerami: {len(total_ltr_elements_with_domains_outside_centromeres)}")

    print(f"Zapisywanie raportu wzbogacenia do: {output_report_file}")
    results_list_domains = []

    with open(output_report_file, 'w') as f_out:
        header_base = [
            "REXdb_Domain_Type",
            "TotalLength_LTRs_In_Centromeres", "TotalLength_LTRs_Outside_Centromeres",
            "Unique_LTR_Count_In_Centromeres", "Unique_LTR_Count_Outside_Centromeres",
            "Prop_Len_Cen", "Prop_Len_NonCen", 
            "Enrich_Factor_Len", "P_Value_Len", "Enriched_Len", "Unique_To_Cen_Len",
            "Prop_Count_Cen", "Prop_Count_NonCen", 
            "Enrich_Factor_Count", "P_Value_Count", "Enriched_Count", "Unique_To_Cen_Count"
        ]
        f_out.write("\t".join(header_base) + "\n")

        for domain_type in sorted(list(all_rexdb_domain_types_found)): # Zmieniono ltr_type na domain_type
            data = domain_data[domain_type] # Użyj domain_type jako klucza
            
            len_cen = data["len_cen"]
            len_noncen = data["len_noncen"]
            other_bp_in_cen_w_domains = max(0, total_ltr_bp_with_domains_in_centromeres - len_cen)
            other_bp_in_noncen_w_domains = max(0, total_ltr_bp_with_domains_outside_centromeres - len_noncen)

            p_val_len, ef_len, enriched_len_str, unique_len_str = float('nan'), 0.0, "Nie", "Nie"
            if total_ltr_bp_with_domains_in_centromeres > 0 or total_ltr_bp_with_domains_outside_centromeres > 0 :
                if len_cen >= 0 and len_noncen >= 0 and other_bp_in_cen_w_domains >= 0 and other_bp_in_noncen_w_domains >= 0:
                    try:
                        table_l = [[int(round(len_cen)), int(round(len_noncen))], 
                                   [int(round(other_bp_in_cen_w_domains)), int(round(other_bp_in_noncen_w_domains))]]
                        if sum(table_l[0]) > 0 and sum(table_l[1]) > 0 and (table_l[0][0] + table_l[1][0]) > 0 and (table_l[0][1] + table_l[1][1]) > 0:
                            _, p_val_len = fisher_exact(table_l, alternative='greater')
                        else: p_val_len = 1.0 
                    except ValueError: pass 

                prop_l_cen = (len_cen / total_ltr_bp_with_domains_in_centromeres) if total_ltr_bp_with_domains_in_centromeres > 0 else 0
                prop_l_noncen = (len_noncen / total_ltr_bp_with_domains_outside_centromeres) if total_ltr_bp_with_domains_outside_centromeres > 0 else 0
                
                if prop_l_noncen > 0: ef_len = prop_l_cen / prop_l_noncen
                elif prop_l_cen > 0: ef_len = float('inf')
                
                if not isinstance(p_val_len, float) or p_val_len >= p_value_threshold_param or len_cen == 0: # Sprawdź, czy p_val_len jest nan
                    enriched_len_str = "Nie"
                elif ef_len > 1 or ef_len == float('inf'): 
                    enriched_len_str = "Tak"
                else:
                    enriched_len_str = "Nie (p<thr, ale EF<=1)"

            unique_len_str = "Tak" if len_cen > 0 and len_noncen == 0 else "Nie"

            count_cen = data["count_cen"] 
            count_noncen = data["count_noncen"]
            other_count_cen = max(0, len(total_ltr_elements_with_domains_in_centromeres) - count_cen)
            other_count_noncen = max(0, len(total_ltr_elements_with_domains_outside_centromeres) - count_noncen)

            p_val_count, ef_count, enriched_count_str, unique_count_str = float('nan'), 0.0, "Nie", "Nie"
            if len(total_ltr_elements_with_domains_in_centromeres) > 0 or len(total_ltr_elements_with_domains_outside_centromeres) > 0:
                if count_cen >= 0 and count_noncen >= 0 and other_count_cen >= 0 and other_count_noncen >= 0:
                    if (count_cen + count_noncen > 0) and (other_count_cen + other_count_noncen > 0) and \
                       (count_cen + other_count_cen > 0) and (count_noncen + other_count_noncen > 0):
                        try:
                            table_c = [[count_cen, count_noncen], [other_count_cen, other_count_noncen]]
                            _, p_val_count = fisher_exact(table_c, alternative='greater')
                        except ValueError: pass
                
                prop_c_cen = (count_cen / len(total_ltr_elements_with_domains_in_centromeres)) if total_ltr_elements_with_domains_in_centromeres else 0
                prop_c_noncen = (count_noncen / len(total_ltr_elements_with_domains_outside_centromeres)) if total_ltr_elements_with_domains_outside_centromeres else 0
                if prop_c_noncen > 0: ef_count = prop_c_cen / prop_c_noncen
                elif prop_c_cen > 0: ef_count = float('inf')

                if not isinstance(p_val_count, float) or p_val_count >= p_value_threshold_param or count_cen == 0:
                    enriched_count_str = "Nie"
                elif ef_count > 1 or ef_count == float('inf'): 
                    enriched_count_str = "Tak"
                else:
                    enriched_count_str = "Nie (p<thr, ale EF<=1)"

            unique_count_str = "Tak" if count_cen > 0 and count_noncen == 0 else "Nie"
            
            row_values = [
                domain_type, len_cen, len_noncen, count_cen, count_noncen,
                f"{prop_l_cen:.4e}", f"{prop_l_noncen:.4e}", f"{ef_len:.2f}", f"{p_val_len:.3e}", enriched_len_str, unique_len_str,
                f"{prop_c_cen:.4e}", f"{prop_c_noncen:.4e}", f"{ef_count:.2f}", f"{p_val_count:.3e}", enriched_count_str, unique_count_str
            ]
            f_out.write("\t".join(map(str,row_values)) + "\n")
            results_list_domains.append({
                "type": domain_type, "len_cen": len_cen, "count_cen": count_cen,
                "ef_len": ef_len, "p_val_len": p_val_len, "enriched_len": enriched_len_str == "Tak",
                "ef_count": ef_count, "p_val_count": p_val_count, "enriched_count": enriched_count_str == "Tak",
            })
        
        f_out.write("\n\n--- PODSUMOWANIE TOP N (Domeny REXdb) ---\n")
        f_out.write(f"(Top {top_n_to_report} typów domen LTR według różnych kryteriów)\n")

        f_out.write(f"\n# Top {top_n_to_report} domen LTR według całkowitej długości LTRów (z tą domeną) w centromerach:\n")
        f_out.write("Domain_Type\tTotalLength_LTRs_In_Centromeres\n")
        for item in sorted(results_list_domains, key=lambda x: x["len_cen"], reverse=True)[:top_n_to_report]:
            if item['len_cen'] > 0:
                f_out.write(f"{item['type']}\t{item['len_cen']}\n")

        f_out.write(f"\n# Top {top_n_to_report} domen LTR według liczby unikalnych LTRów (z tą domeną) w centromerach:\n")
        f_out.write("Domain_Type\tUnique_LTR_Count_In_Centromeres\n")
        for item in sorted(results_list_domains, key=lambda x: x["count_cen"], reverse=True)[:top_n_to_report]:
            if item['count_cen'] > 0:
                f_out.write(f"{item['type']}\t{item['count_cen']}\n")
        
        f_out.write(f"\n# Top {top_n_to_report} domen LTR według współczynnika wzbogacenia (długość, p < {p_value_threshold_param}):\n")
        f_out.write("Domain_Type\tEnrichment_Factor_Length\tP_Value_Length\n")
        enriched_len_items_domain = [item for item in results_list_domains if item['enriched_len']]
        for item in sorted(enriched_len_items_domain, key=lambda x: x["ef_len"], reverse=True)[:top_n_to_report]:
            f_out.write(f"{item['type']}\t{item['ef_len']:.2f}\t{item['p_val_len']:.3e}\n")

        f_out.write(f"\n# Top {top_n_to_report} domen LTR według współczynnika wzbogacenia (liczba LTRów, p < {p_value_threshold_param}):\n")
        f_out.write("Domain_Type\tEnrichment_Factor_Count\tP_Value_Count\n")
        enriched_count_items_domain = [item for item in results_list_domains if item['enriched_count']]
        for item in sorted(enriched_count_items_domain, key=lambda x: x["ef_count"], reverse=True)[:top_n_to_report]:
            f_out.write(f"{item['type']}\t{item['ef_count']:.2f}\t{item['p_val_count']:.3e}\n")

    print(f"Analiza zakończona. Raport zapisano w: {output_report_file}")

if __name__ == '__main__':
    # Sprawdzenie czy translate_seq jest dostępne, jeśli Etap 1 ma być wykonany
    # (ta logika jest bardziej w run_hmmsearch_on_ltrs)
    if not os.path.exists(os.path.join(SCRIPT_PATH_BASE, "translate_seq.py")):
         # To ostrzeżenie jest ważne, jeśli użytkownik nie ma jeszcze wyników Etapu 1
         print(f"OSTRZEŻENIE: Nie znaleziono pliku translate_seq.py w katalogu {SCRIPT_PATH_BASE}")
         print("             Jeśli Etap 1 (analiza HMMER) nie został jeszcze pomyślnie wykonany,")
         print("             ten plik będzie potrzebny. Jeśli Etap 1 jest pomijany, to ostrzeżenie można zignorować.")
    main()
