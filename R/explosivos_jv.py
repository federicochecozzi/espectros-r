import matplotlib.pyplot as plt
from matplotlib import interactive#necesario si estoy usando IDLE
interactive(True)
import os
from fnmatch import fnmatch
import ntpath
import numpy as np

# Agregado Juan:
a= os.getcwd()
directorio = a + "\\datos"
archivos = os.listdir(directorio)       
# Modificado Juan:
# archivos =          ['PE1.csv', 'PE2.csv', 'PE3.csv', 'Pu2.csv', 'Pu3.csv', 'RD1.csv', 'RD2.csv', 'RD3.csv', 'TN1.csv', 'TN2.csv', 'TN3.csv']
# archivos = obtener_archivos(directorio)
espectros_validos = [ [3,4,5],   [3,4,5],   [2,3,4],   [2,3,4],   [2,3,4],   [2,3,4],   [2,3,4],   [2,3,4],   [2,3,4],    [2,3,4],    [2,3,4]]

integrales = [[471,474.5], [488,497], [498,506], [509,515]]
picos = [472.6, 492.99, 500.7, 512.24]


def normalizar_espectros(long_ondas, espectro, longitud_de_onda):
    indice = min(range(len(long_ondas)), key=lambda i: abs(long_ondas[i]-longitud_de_onda))
    val = espectro[indice]
    espectro_normalizado = []
    for i in espectro:
        espectro_normalizado.append(i/val)
    return espectro_normalizado

def recortar_espectro(long_ondas, espectro, long_onda_a, long_onda_b):
    indice_a = min(range(len(long_ondas)), key=lambda i: abs(long_ondas[i]-long_onda_a))
    indice_b = min(range(len(long_ondas)), key=lambda i: abs(long_ondas[i]-long_onda_b))
    return long_ondas[indice_a:indice_b], espectro[indice_a:indice_b]


def obtener_archivos(directorio):
    archivos = [] #Lista de archivos
    for path, subdirs, files in os.walk(directorio): #Busca todos los archivos que coincidan con la extension el el directorio
        for name in files:
            if fnmatch(name, '.csv'):
                archivos.append(os.path.join(path, name))
    return archivos


def leer_espectro_csv( filename):
    spectra = []
    wave_length = []
    f = open(filename, "r")
    f.readline()
    f.readline()
    line = f.readline()
    while line != '':
        line = line.replace(',', '.')
        splitted_line = line.split(';')
        while(len(spectra) < len(splitted_line) - 2):
            spectra.append([])
        wave_length.append(float(splitted_line[0]))
        i = 0
        for intensity in splitted_line[1:-1]:
            spectra[i].append(float(intensity))
            i += 1
        line = f.readline()
    f.close()
    return wave_length, spectra

def calcular_integral(long_ondas, espectro, long_onda_a, long_onda_b):
    indice_a = min(range(len(long_ondas)), key=lambda i: abs(long_ondas[i]-long_onda_a))
    indice_b = min(range(len(long_ondas)), key=lambda i: abs(long_ondas[i]-long_onda_b))
    integral = 0.0
    for i in range(indice_a, indice_b):
        integral += espectro[i] * (long_ondas[i+1] - long_ondas[i])
    return integral

def valor(long_ondas, espectro, long_onda):
    indice = min(range(len(long_ondas)), key=lambda i: abs(long_ondas[i]-long_onda))
    return espectro[indice]

all_wave_lengths = []
all_spectra = [] 
"""
[
    [[espectro],[espectro],[espectro],...], #tnt1.csv

    [[espectro],[espectro],[espectro],...], #tnt2.csv

    [[espectro],[espectro],[espectro],...], #tnt3.csv

    ...
]
"""
               
for archivo in archivos:
    # Modificado Juan:
    archivo =  directorio + "\\" + archivo
    wave_length, spectra = leer_espectro_csv(archivo)
    all_spectra.append(spectra)
    all_wave_lengths.append(wave_length)

n_espectro = 0

# Modificado Juan:
f_out = open("resultados1.csv", 'w+')
f_out.write(" ,")
for integral in integrales:
    f_out.write("Integral " + str(integral[0]) + " a " + str(integral[1]) + ",")
for pico in picos:
    f_out.write("Pico " + str(pico) + ",")

f_out.write('\n')

for spectra in all_spectra:

    for i in espectros_validos[n_espectro]:

        spec_norm = normalizar_espectros(all_wave_lengths[0], spectra[i], 481.8)

        wave_length_cut, spec_cut = recortar_espectro(all_wave_lengths[0], spec_norm, 465, 520)
      # Modificado Juan:
        f_out.write(archivos[n_espectro][:-4] + "-" + str(i) + ' ,')

        for integral in integrales:
            val_integral = calcular_integral(wave_length_cut, spec_cut, integral[0], integral[1])
            print(archivos[n_espectro] + " - " + str(i) + " - Integral de " + str(integral[0]) + " a " + str(integral[1]) + " = " + str(val_integral))
            f_out.write('%.4f,' % val_integral)

        for pico in picos:
            val = valor(wave_length_cut, spec_cut, pico)
            f_out.write('%.4f,' % val)
        plt.plot(wave_length_cut, spec_cut)

        f_out.write('\n')

    print("")

    n_espectro += 1

f_out.close()
