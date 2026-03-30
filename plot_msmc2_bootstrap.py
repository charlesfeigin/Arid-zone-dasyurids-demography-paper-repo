import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['svg.fonttype'] = 'none'

mu = 5.95e-9
gen = 1

def load_bootstraps(mainfile, bootprefix):
    dfs = []

    # main
    dfs.append(pd.read_csv(mainfile, sep=r"\s+"))

    # bootstraps 1 to 100
    for i in range(1, 101):
        fn = f"{bootprefix}.bootstrap_{i}.final.txt"
        dfs.append(pd.read_csv(fn, sep=r"\s+"))

    return dfs

Smurina = load_bootstraps(
    "slendertailed/whole_genome_msmc2_slendertailed_i50_p1-2_25-1_1-2_1-3.final.txt",
    "slendertailed/slendertailed_i50_p1-2_25-1_1-2_1-3"
)

Nridei = load_bootstraps(
    "wongai/whole_genome_msmc2_wongai_i50_p1-2_25-1_1-2_1-3.final.txt",
    "wongai/wongai_i50_p1-2_25-1_1-2_1-3"
)

for df in Smurina[1:]:
    plt.step(df["left_time_boundary"]/mu*gen, (1/df["lambda"])/(2*mu),
             color="violet", linewidth=0.75, alpha=0.1)

plt.step(Smurina[0]["left_time_boundary"]/mu*gen, (1/Smurina[0]["lambda"])/(2*mu),
         color="darkviolet", linewidth=1.75)

for df in Nridei[1:]:
    plt.step(df["left_time_boundary"]/mu*gen, (1/df["lambda"])/(2*mu),
             color="green", linewidth=0.75, alpha=0.1)

plt.step(Nridei[0]["left_time_boundary"]/mu*gen, (1/Nridei[0]["lambda"])/(2*mu),
         color="darkgreen", linewidth=1.75)

plt.xlim(0, 250000)
plt.ylim(0, 1000000)

# Area of higher confidence (1000g/y bp to 100kg/ky bp)
plt.axvline(x=1000, color="black", alpha=0.75, linewidth=1, linestyle='--')
plt.axvline(x=100000, color="black", alpha=0.75, linewidth=1, linestyle='--')


# Beginning of humid phase in murry-darling basin
plt.axvline(x=40000, color="blue", alpha=1, linewidth=2.0)

# Beginning of increased fluvial activity in murray-darling basin
plt.axvline(x=35000, color="darkblue", alpha=1, linewidth=2.0)


# Decreased flows in murray
plt.axvspan(xmin=7000, xmax=14000, alpha=0.5, color='pink')

plt.tick_params(axis="both", pad=16)

ax = plt.gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(0.5)
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(0.5)

plt.savefig("plot.svg", format="svg")
