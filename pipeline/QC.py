import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def qc_plots(df, umbral=20, trimmed=True):
    
    sns.set(style="whitegrid")
    figs = []

    #Se toman los datos de los dataframes
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

    df = df.copy()
    df["read_ID"] = (df["ID"].astype(str).str.replace(".ab1", "", regex=False).str.split("_").str[:2].str.join("_"))
    df["Sample"] = df["read_ID"].str.split("_").str[0]
    df = df.sort_values(["Sample", "read_ID"])

    #Longitud de las secuencias.
    fig1, ax1 = plt.subplots(figsize=(10, 5))
    sns.barplot(data=df, x="read_ID", y=len_col, ax=ax1)
    ax1.axhspan(0, 400, color="red", alpha=0.15)
    ax1.axhspan(400, 650, color="orange", alpha=0.15)
    ax1.axhspan(650, df[len_col].max() * 1.1, color="green", alpha=0.15)
    ax1.set_title("Longitud de secuencia")
    ax1.tick_params(axis="x", rotation=90, labelsize=6)
    figs.append(fig1)
    
    #Boxplot con calidades de las secuencias.
    cal_data = []
    for _, row in df.iterrows():
        for q in row[cal_col]:
            cal_data.append({"Sample": row["Sample"], "read_ID": row["read_ID"], "Quality": q})

    cal_df = pd.DataFrame(cal_data)

    fig2, ax2 = plt.subplots(figsize=(10, 5))
    sns.boxplot(data=cal_df, x="read_ID", y="Quality", ax=ax2)
    ax2.axhspan(0, 20, color="red", alpha=0.15)
    ax2.axhspan(20, 30, color="orange", alpha=0.15)
    ax2.axhspan(30, cal_df["Quality"].max() * 1.1, color="green", alpha=0.15)
    ax2.set_title("Distribución de calidades")
    ax2.tick_params(axis="x", rotation=90, labelsize=6)
    figs.append(fig2)

    #Media y mediana de las secuencias.
    fig3, ax3 = plt.subplots(figsize=(6, 5))

#Scatterplot sin leyenda
    sns.scatterplot(
        data=df,
        x=mean_col, y=median_col,
        s=70, ax=ax3,
        color="steelblue")

    #Líneas de umbral
    ax3.axvline(x=umbral, color="purple", linestyle="--", linewidth=1.25)
    ax3.axhline(y=umbral, color="purple", linestyle="--", linewidth=1.25)

    #Zona de buena calidad.
    ax3.fill_betweenx(
    y=[umbral, ax3.get_ylim()[1]],
    x1=umbral, x2=ax3.get_xlim()[1],
    color="green",alpha=0.1)

# Se etiquetan las muestras problemáticas.
    bad = df[ (df[mean_col] < umbral) | (df[median_col] < umbral)]

    for _, row in bad.iterrows():
        ax3.text(
            row[mean_col],
            row[median_col],
            row["Sample"],
            fontsize=7, alpha=0.8,
            ha="right", va="bottom")

    #Etiquetas y título
    ax3.set_title("Calidad media vs mediana")
    ax3.set_xlabel("Calidad media")
    ax3.set_ylabel("Calidad mediana")

    figs.append(fig3)

    return figs
