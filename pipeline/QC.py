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
    df["Sample"] = df["ID"].astype(str)

    #Longitud de las secuencias.
    fig1, ax1 = plt.subplots(figsize=(10, 4))
    sns.barplot(data=df, x="Sample", y=len_col, ax=ax1)
    ax1.axhspan(0, 400, color="red", alpha=0.15)
    ax1.axhspan(400, 650, color="orange", alpha=0.15)
    ax1.axhspan(650, df[len_col].max() * 1.1, color="green", alpha=0.15)
    ax1.set_title("Longitud de secuencia")
    ax1.tick_params(axis="x", rotation=45)
    figs.append(fig1)
    
    #Boxplot con calidades de las secuencias.
    cal_data = []
    for _, row in df.iterrows():
        for q in row[cal_col]:
            cal_data.append({"Sample": row["Sample"], "Quality": q})

    cal_df = pd.DataFrame(cal_data)

    fig2, ax2 = plt.subplots(figsize=(12, 5))
    sns.boxplot(data=qual_df, x="Sample", y="Quality", ax=ax2)
    ax2.axhspan(0, 20, color="red", alpha=0.15)
    ax2.axhspan(20, 30, color="orange", alpha=0.15)
    ax2.axhspan(30, cal_df["Quality"].max() * 1.1, color="green", alpha=0.15)
    ax2.set_title("Distribución de calidades")
    ax2.tick_params(axis="x", rotation=45)
    figs.append(fig2)

    #Media y mediana de las secuencias.
    fig3, ax3 = plt.subplots(figsize=(6, 5))
    sns.scatterplot(data=df, 
                    x=mean_col, y=median_col,
                    hue="Sample", s=80,ax=ax3)
    ax3.axvline(x=umbral,
                color="red",
                linestyle="--", linewidth=1.5,
                label=f"Umbral = {umbral}")
    ax3.axhline(y=umbral,
                color="red",
                linestyle="--", linewidth=1.5)
    ax3.fill_betweenx(y=[threshold, ax3.get_ylim()[1]], 
                      x1=threshold, x2=ax3.get_xlim()[1], 
                      color="green", alpha=0.1)

    ax3.set_title("Calidad media vs mediana")
    ax3.set_xlabel("Calidad media")
    ax3.set_ylabel("Calidad mediana")

    handles, labels = ax3.get_legend_handles_labels()
    ax3.legend(handles[:len(df["Sample"].unique()) + 1], labels[:len(df["Sample"].unique()) + 1])

    figs.append(fig3)

    return figs
