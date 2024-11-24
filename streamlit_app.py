import streamlit as st
from Bio.Seq import Seq
from Bio import SeqIO

# Título de la aplicación
st.title("Contador de Proteínas en Secuencia de ADN")

# Subtítulo
st.header("Ingresa una secuencia de ADN para contar las proteínas y mostrar sus nombres")

# Entrada de la secuencia de ADN
adn_sequence = st.text_area("Secuencia de ADN:")

# Función para traducir la secuencia de ADN a proteínas
def translate_to_proteins(seq):
    # Asegurarse de que la longitud sea múltiplo de 3
    if len(seq) % 3 != 0:
        st.warning("La longitud de la secuencia no es múltiplo de 3, lo que puede generar proteínas incompletas.")
    
    # Traducir la secuencia de ADN en fragmentos de 3 nucleótidos
    proteins = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if len(codon) == 3:
            # Traducción del codón en proteína
            protein = str(Seq(codon).translate())
            proteins.append(protein)
    
    return proteins

# Función para contar las proteínas y mostrar sus nombres
def count_and_show_proteins(seq):
    proteins = translate_to_proteins(seq)
    if not proteins:
        st.write("No se encontraron proteínas en la secuencia.")
        return
    
    # Mostrar el número de proteínas
    st.write(f"La secuencia contiene {len(proteins)} proteína(s).")
    
    # Mostrar los nombres (secuencias) de las proteínas
    for i, protein in enumerate(proteins):
        st.write(f"Proteína {i+1}: {protein}")

# Si el usuario ha ingresado una secuencia de ADN, realizar el análisis
if adn_sequence:
    count_and_show_proteins(adn_sequence)
