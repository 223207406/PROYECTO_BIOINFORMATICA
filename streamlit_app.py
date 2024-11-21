import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

st.title("¿Cuántas proteínas contiene y cuáles son? JUAN RAMONSITO DE CLN")
st.write(
    "Let's start building! For help and inspiration, head over to [docs.streamlit.io](https://docs.streamlit.io/)."
)

# Diccionario para la traducción del código genético
codon_to_amino_acid = {
    'AUG': 'Met', 'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu', 'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met', 'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser', 'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr', 'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop', 'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys', 'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp', 'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg', 'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}

# Función para transcribir ADN a ARN
def transcribir_adn_a_arn(adn):
    return adn.replace('T', 'U')

# Función para traducir ARN a proteína
def traducir_arn_a_proteina(arn):
    proteina = []
    for i in range(0, len(arn), 3):
        codon = arn[i:i+3]
        if len(codon) == 3:
            amino_acido = codon_to_amino_acid.get(codon, 'Stop')
            if amino_acido == 'Stop':
                break
            proteina.append(amino_acido)
    return proteina

# Función para contar la cantidad de tipos de proteínas
def contar_proteinas(proteina):
    return dict(Counter(proteina))

# Función para graficar las proteínas
def graficar_proteinas(conteo_proteinas):
    # Datos para la gráfica
    proteinas = list(conteo_proteinas.keys())
    cantidades = list(conteo_proteinas.values())

    # Crear el gráfico
    plt.figure(figsize=(10, 6))
    plt.bar(proteinas, cantidades, color='skyblue')

    # Añadir título y etiquetas
    plt.title('Frecuencia de Proteínas en la Secuencia de ADN')
    plt.xlabel('Proteínas')
    plt.ylabel('Frecuencia')
    plt.xticks(rotation=45, ha="right")

    # Mostrar el gráfico
    plt.tight_layout()
    plt.show()

# Función principal
def analizar_adn_y_graficar(adn):
    # Transcripción de ADN a ARN
    arn = transcribir_adn_a_arn(adn)
    
    # Traducción de ARN a proteína
    proteina = traducir_arn_a_proteina(arn)
    
    # Contar los tipos de proteínas
    conteo_proteinas = contar_proteinas(proteina)
    
    # Graficar los resultados
    graficar_proteinas(conteo_proteinas)

secuencia_adn = "ATGCGTACGTTAGC" 

# Analizar la secuencia de ADN y graficar
analizar_adn_y_graficar(secuencia_adn)
