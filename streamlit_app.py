import streamlit as st

st.title("Proyecto Bioinformática.")
st.write(
    "Let's start building! For help and inspiration, head over to [docs.streamlit.io](https://docs.streamlit.io/)."
)
import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight
from Bio import pairwise2

# Título de la aplicación
st.title("Comparador de Secuencias de ADN y Proteínas: E. Coli vs Salmonella")

# Subtítulo
st.header("Ingresa las secuencias de ADN de E. Coli y Salmonella para compararlas")

# Entrada de secuencias de ADN
ecoli_sequence = st.text_area("Secuencia de ADN de E. Coli:")
salmonella_sequence = st.text_area("Secuencia de ADN de Salmonella:")

# Función para analizar las secuencias
def analyze_sequences(seq1, seq2):
    # Asegurarnos de que las secuencias sean válidas
    if len(seq1) == 0 or len(seq2) == 0:
        st.warning("Por favor ingresa ambas secuencias de ADN.")
        return

    # Longitud de las secuencias
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)
    
    # Contar las bases A, T, C, G para cada secuencia
    bases_ecoli = {base: seq1.count(base) for base in "ATCG"}
    bases_salmonella = {base: seq2.count(base) for base in "ATCG"}

    # Comparar las secuencias usando alineamiento global
    alignments = pairwise2.align.globalxx(seq1, seq2)
    best_alignment = alignments[0]
    match_score = best_alignment[2]  # Puntaje de coincidencia
    
    # Mostrar resultados
    st.subheader("Resultados de la Comparación:")
    st.write(f"Longitud de la secuencia E. Coli: {len_seq1} nucleótidos")
    st.write(f"Longitud de la secuencia Salmonella: {len_seq2} nucleótidos")
    
    st.write(f"Porcentaje de coincidencia entre E. Coli y Salmonella: {match_score / max(len_seq1, len_seq2) * 100:.2f}%")

    st.write("Bases en la secuencia E. Coli:")
    st.write(bases_ecoli)
    st.write("Bases en la secuencia Salmonella:")
    st.write(bases_salmonella)

    st.write(f"\nEl puntaje de alineamiento global es: {match_score}")

    # Mostrar las proteínas codificadas por cada secuencia
    ecoli_proteins = translate_to_proteins(seq1)
    salmonella_proteins = translate_to_proteins(seq2)
    
    st.subheader("Proteínas Codificadas:")
    st.write(f"E. Coli tiene {len(ecoli_proteins)} proteínas:")
    for i, protein in enumerate(ecoli_proteins):
        st.write(f"Proteína {i+1}: {protein}")
    
    st.write(f"Salmonella tiene {len(salmonella_proteins)} proteínas:")
    for i, protein in enumerate(salmonella_proteins):
        st.write(f"Proteína {i+1}: {protein}")

# Función para traducir una secuencia de ADN a proteínas
def translate_to_proteins(seq):
    # Asegurarnos de que la longitud de la secuencia sea múltiplo de 3 para traducción
    if len(seq) % 3 != 0:
        st.warning("La secuencia no tiene una longitud múltiplo de 3, lo que puede generar proteínas incompletas.")

    # Traducir la secuencia en trozos de 3 nucleótidos (códones)
    proteins = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if len(codon) == 3:
            protein = str(Seq(codon).translate())
            proteins.append(protein)
    
    return proteins

# Si el usuario ha ingresado ambas secuencias, analizamos
if ecoli_sequence and salmonella_sequence:
    analyze_sequences(ecoli_sequence, salmonella_sequence)
