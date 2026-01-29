import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def qc_plots(df, umbral=20, trimmed=True):

    sns.set(style="whitegrid")
    figs = []

    #Se toman las columnas según estado.
    if trimmed:
        seq_col = "Secuencia cortada"
        cal_col = "Puntuación cortada"
        len_col = "Longitud (pb)"
        mean_col = "Calidad Media"
        median_col = "Mediana calidad"
        pct_col = "% Bases alta calidad"
        amb_col = "Bases ambiguas (N)"
    else:
        seq_col = "Secuencia Fasta"
        cal_col = "Puntuación calidad"
        len_col = "Longitud (pb)"
        mean_col = "Calidad media"
        median_col = "Mediana calidad"
        pct_col = "% Alta calidad"
        amb_col = "Bases ambiguas(N)"

    #Se preparan los identificadores.
    df_qc = df.copy()

    #Se identifican las muestras y sus primers.
    df_qc["read_ID"] = (df["ID"].astype(str).str.replace(".ab1", "", regex=False).str.split("_").str[:2].str.join("_"))
    df_qc["Sample"] = df_qc["read_ID"].str.split("_").str[0]
    df_qc = df_qc.sort_values(["Sample", "read_ID"]).reset_index(drop=True)

    #Se identifica la posición en que cambia la muestra.
    cambio = df_qc["Sample"].ne(df_qc["Sample"].shift())
    separador = [i - 0.5 for i in df.index[cambio][1:]]

    #Histograma de Longitud de secuencias.
    fig1, ax1 = plt.subplots(figsize=(12, 5))

    sns.barplot(data=df_qc,x="read_ID",y=len_col,ax=ax1)

    ax1.axhspan(0, 400, color="red", alpha=0.15)
    ax1.axhspan(400, 650, color="orange", alpha=0.15)
    ax1.axhspan(650, df[len_col].max() * 1.1, color="green", alpha=0.15)

    for pos in separador:
        ax1.axvline(pos, color="black", linestyle="-", linewidth=0.6, alpha=0.5)

    ax1.set_title("Longitud de secuencia")
    ax1.tick_params(axis="x", rotation=90, labelsize=5)

    figs.append(fig1)

    #Boxplot con las calidades de las secuencias.
    cal_data = []
    for _, row in df_qc.iterrows():
        for q in row[cal_col]:
            cal_data.append({
                "read_ID": row["read_ID"],
                "Sample": row["Sample"],
                "Quality": q})

    cal_df = pd.DataFrame(cal_data)

    fig2, ax2 = plt.subplots(figsize=(12, 5))

    sns.boxplot(data=cal_df, x="read_ID",y="Quality",ax=ax2)
    ax2.axhspan(0, 20, color="red", alpha=0.15)
    ax2.axhspan(20, 30, color="orange", alpha=0.15)
    ax2.axhspan(30, cal_df["Quality"].max() * 1.1, color="green", alpha=0.15)

    for pos in separador:
        ax2.axvline(pos, color="black", linestyle="-", linewidth=0.6, alpha=0.5)

    ax2.set_title("Distribución de calidades")
    ax2.tick_params(axis="x", rotation=90, labelsize=5)

    figs.append(fig2)

    #Media vs mediana de las calidades.
    fig3, ax3 = plt.subplots(figsize=(6, 5))

    sns.scatterplot(data=df_qc,
        x=mean_col, y=median_col,
        ax=ax3, s=70, color="steelblue")

    ax3.axvline(x=umbral, color="red", linestyle="--", linewidth=1.5)
    ax3.axhline(y=umbral, color="red", linestyle="--", linewidth=1.5)

    ax3.fill_betweenx(
        y=[umbral, ax3.get_ylim()[1]],
        x1=umbral, x2=ax3.get_xlim()[1],
        color="green", alpha=0.1)

    ax3.set_title("Calidad media vs mediana")
    ax3.set_xlabel("Calidad media")
    ax3.set_ylabel("Calidad mediana")

    figs.append(fig3)

    return figs
