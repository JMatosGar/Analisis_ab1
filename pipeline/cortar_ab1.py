import statistics

def cortar_ab1(df, umbral=20, min_bases=5, output="_trimmed"):
    resultados = []

    for _, row in df.iterrows():
        secuencia = row["Secuencia Fasta"]
        calidad = row["Puntuación calidad"]
        sec_len = len(secuencia)

        if sec_len == 0 or not calidad:
            continue

      #Se identifican los puntos del inicio y fin de la secuencia donde la calidad es suficiente.
        inicio = 0
        for i in range(sec_len - min_bases + 1):
            if all(q >= umbral for q in calidad[i:i + min_bases]):
                inicio = i
                break

        fin = sec_len
        for i in range(sec_len, min_bases - 1, -1):
            if all(q >= umbral for q in calidad[i - min_bases:i]):
                fin = i
                break

        # ---- Apply trimming ----
        trimmed_sec = secuencia[inicio:fin]
        trimmed_cal = calidad[inicio:fin]

        trimmed_len = len(trimmed_sec)
        bases_trimmed = sec_len - trimmed_len

        if trimmed_len == 0:
            continue

        #Se recalculan los parámetros.
        calidad_media = statistics.mean(trimmed_cal)
        calidad_mediana = statistics.median(trimmed_cal)

        calidad_alta = sum(1 for q in trimmed_cal if q >= umbral)
        pct_calidad_alta = (calidad_alta / trimmed_len) * 100

        num_amb = trimmed_sec.upper().count("N")

        larga_calidad = 0
        actual = 0
        for q in trimmed_cal:
            if q >= umbral:
                actual += 1
                larga_calidad = max(larga_calidad, actual)
            else:
                actual = 0

        resultados.append({
            "ID": row["ID"],
            "Secuencia cortada": trimmed_sec,
            "Puntuación cortada": trimmed_cal,
            "Longitud (pb)": trimmed_len,
            "Bases cortadas": bases_trimmed,
            "Inicio": inicio,
            "Fin": fin,
            "Calidad Media": calidad_media,
            "Mediana calidad": calidad_mediana,
            "% Bases alta calidad": pct_calidad_alta,
            "Bases ambiguas (N)": num_amb,
            "Tramo alta calidad": larga_calidad})

    return trimmed_df
