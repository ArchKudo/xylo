- Abstract (500)
- Acknowledgements (-)
- Table of Contents (-)
- List of Figures
- Introduction (2000)
  - Why Xylitol
  - Why Xylose reductase
  - Why candida troicalis
  - Why cattle rumen
  - Background
    - Gene
    - Chromosome
    - Alleles
    - Haplotypes
    - Cattle Rumen
    - Stereoisomers
    - Xylitol
    - Genome/Metagenome
  - Motivation
  - Aims
  - Objectives

- Who is the customer?
  - LTS customer: People consuming artificial sweeteners
  - Immediate customer: Wet lab tech, who can culture the provided haplotype
- What’s the problem you’re trying to solve?
  - Find if an organism exists which produces XR in the microbial community of a cattle's rumen
- What’s the solution (and you need to explain to the customer, not to yourself).
  - Use metagenomics
    - Get known XR gene variants
    - Align it to the metagenome reads of the rumen
    - Find variants
    - Construct the haplotype
- Would they reasonably adopt this solution, because you’re asking for a behaviour change?
  - Known procedure exist to culture only xylitol and not rabinol using the dark arts, allowing for easier onboarding
- What’s the Total Addressable Market (TAM), and is it big enough to be worth doing?
  - Can expand to research/industry like as it provides reduced cost of production

Introduction

~~Xylitol is one of the naturally occurring pentitols (five-carbon sugar alcohol) with a molecular formula of C5H12O5. It is a white crystalline sugar, commercially used as an artificial sweetener in the food and pharmaceutical industries. Xylitol has a lesser calorific value (2.4 cal/g) compared to sucrose (4.0 cal/g) but has a relative sweetness almost equal to sucrose (Chen et al., 2010; Tiefenbacher, 2017). The xylitol metabolism is independent of insulin; hence it has been used as a safe sucrose substitute for patients with diabetes. In the human gastrointestinal tract, 50–75% ingested xylitol is not absorbed (Rehman et al., 2013). Owing to its anti-cariogenic properties, it has been used in the production of chewing gums and toothpaste. Many studies show the role of xylitol to reduce the incidence of respiratory and middle ear infections. Xylitol increases the absorption of calcium thereby helps in combating osteoporosis (Mussatto, 2012). In 1891, a German and French chemist concurrently discovered xylitol (Rehman et al., 2016). Several years later, xylitol crystals were purified successfully, and it was used widely as an alternative to sugar during World War II, to meet the severe shortage faced during the war (Rehman et al., 2013). Xylitol was recognized as a potent sweetener after the discovery of its property of not elevating blood sugar levels. Countries like Germany, Japan, and Italy started to use xylitol as a part of the diabetic diet. Later in the 1970s, the dental benefits of xylitol came into light after the “Turku Sugar Studies” (Janakiram et al., 2017). Since 1975, various xylitol substituted products were introduced globally. ylitol is naturally present in little quantity in diverse vegetables and fruits. Some of the fruits containing xylitol are strawberry, raspberry, banana, yellow plums, and vegetables containing xylitol include cauliflower, spinach, carrot, onion, white mushroom, eggplant, lettuce, and pumpkin. It is also found in hardwood trees like birch and beechwood and in husks and stalks of the plants (Chen et al., 2010). Humans and animals produce small quantities of xylitol as well, during the metabolism of glucose. On average, an adult human produces 5 to 15 g of xylitol per day (Rehman et al., 2016). Since the amount of xylitol present in these natural sources is low, its extraction from these sources is inefficient. Currently, the chemical method of xylitol production is being used commercially to meet the xylitol demand. In this chemical process, xylose extracted from the lignocellulosic biomass undergoes catalytic hydrogenation to produce xylitol. Alternatively, to reduce the high production costs, biotechnological methods of xylitol production from lignocellulosic biomass can be employed (Mathew et al., 2018). There is an increasing interest in commercializing biotechnological methods for the sustainable production of xylitol.~

Xylitol is a sugar alcohol which can be used as an alternative to sugar in food[Umai, 2022] as well as in pharmaceutical products like chewing gum used for treatment of ear infections and dentistry due to its anticariogenic properties.[2, find more uses] and is also known to increase calcium absorption. Large scale production of xylitol is currently achieved using industrial chemical processing. The aim of this project is to find biological enzymes which can sustainably help in reducing the energy costs of production

~~In recent years, advancements in computational biology and data science have opened new avenues for exploring and understanding the biological world. One of the promising areas of research is metagenomics, which involves the study of DNA sequences from diverse organisms co-existing within specific environments. Metagenomics allows researchers to delve into complex microbial communities that play crucial roles in various ecological niches, including the breakdown of biomass, antimicrobial defense, and adaptation to extreme conditions. Within these communities, enzymes, the catalysts of biochemical reactions, play an essential role in fulfilling specific biological functions.

This research proposal seeks to compare and contrast two distinct approaches for studying enzymes from the aldo/keto reductase (AKR2) family, which are commercially important and have significant applications in diverse industries. The first approach involves the computational recovery of haplotype enzymes from metagenomic DNA sequences using the Hansel and Gretel algorithm. The second approach utilizes protein embeddings generated by large language models, such as the ProtT5 Model, to represent the structural and functional properties of AKR2 enzymes.~~

Background

To better understand the proposed research, it is essential to establish a foundational knowledge of key biological concepts and terminologies related to genetics, metagenomics, and protein embeddings.

~~Gene, Chromosome, Alleles, and Haplotypes:
Genes are segments of DNA that contain instructions for the synthesis of specific proteins or functional RNA molecules. Genes are organized into chromosomes, and each chromosome carries numerous genes. In diploid organisms, including humans, there are two copies of each chromosome, and consequently, two copies of each gene, referred to as alleles. Haplotypes are combinations of alleles that are inherited together from one parent, and they play a crucial role in understanding the genetic variation within a population.

Cattle Rumen:
The cattle rumen is the first compartment of the stomach in cattle and other ruminant animals. It hosts a complex microbial ecosystem that efficiently breaks down plant materials through the fermentation process. This environment is an example of a metagenomic niche rich in enzymes with diverse functions.

Stereoisomers and Xylitol:
Stereoisomers are molecules that have the same chemical formula and bonding pattern but differ in their three-dimensional arrangement. Xylitol is a specific example of a stereoisomer and is a sugar alcohol used as a sugar substitute. Enzymes play a critical role in the synthesis and degradation of such compounds.

Genome/Metagenome:
A genome refers to the complete set of genetic material within an organism, while a metagenome represents the collective genetic material of all the organisms present in a specific environmental sample. Metagenomics involves the sequencing and analysis of metagenomes to gain insights into the microbial community and its functional potential.~~

Motivation

~~The motivation behind this research stems from the importance of understanding the functional diversity of enzymes in the AKR2 family, which has significant relevance in various industrial applications. By investigating and comparing the computational recovery of haplotype enzymes and the utilization of protein embeddings, we aim to shed light on the strengths and limitations of each approach.

Metagenomics offers a powerful means to explore the vast reservoir of genetic diversity within microbial communities. However, the analysis of metagenomic data is inherently challenging due to the complexity of the samples and the co-existence of multiple organisms. The Hansel and Gretel algorithm addresses this challenge by computationally recovering haplotype enzymes, providing valuable insights into the genetic variations within enzyme families.

On the other hand, protein embeddings offer an alternative method to represent and characterize enzymes in a machine-friendly format. These embeddings provide compact, vector-based representations of protein sequences, enabling efficient processing and analysis. By leveraging large language models, such as the ProtT5 Model, we can extract valuable structural and functional information from AKR2 enzymes.~~

Aims

~~The primary aim of this research project is to compare and contrast the haplotype-based computational approach with the protein embedding-based approach in characterizing AKR2 enzymes recovered from metagenomic DNA sequences.

To achieve this, the specific aims of the study are as follows:

    Utilize the Hansel and Gretel algorithm to computationally recover haplotype enzymes from metagenomic DNA sequences sampled from the cattle rumen and other relevant environments.

    Leverage the ProtT5 Model, a large language model, to generate protein embeddings that capture the structural and functional properties of AKR2 enzymes.

    Perform a comprehensive comparative analysis of the results obtained from the two approaches, focusing on the ability of each method to identify and analyze specific variations within the AKR2 enzyme family.

    Evaluate the effectiveness of protein embeddings in representing and characterizing AKR2 enzymes compared to the computational recovery of haplotype enzymes.

    Gain insights into the strengths, limitations, and potential applications of both methods for studying commercially important enzymes and their relevance in various industries.

In conclusion, this research proposal aims to contribute to the field of metagenomics and protein bioinformatics by providing a comprehensive analysis of two distinct approaches for studying the AKR2 enzyme family. The findings of this study have the potential to impact multiple industries that rely on enzymes for various applications, ranging from biotechnology to agriculture and beyond. Moreover, understanding the strengths and limitations of each approach will pave the way for future advancements in enzyme research and facilitate the design of novel biocatalysts with enhanced efficiency and specificity.~~


- Literature Review (2000)
  - The literature review is all about the related knowledge that you are building on.  Similar products and related research are usual. Remember to use your own words and to show relevance to your project aim. The literature review will refer extensively to the bibliography.  Harvard (author-date) and IEEE reference styles are usual in Computer Science, but the only real rule is that you should use a consistent style.


    
- Report (8000)
  - Reporting on the project will normally require more than one chapter. A development project is likely to have chapters addressing requirements, design, implementation, testing and packaging if a plan-driven method is used.  If an agile approach is taking, you might have a chapter for each sprint or iteration.
  - Hungate
  - Descriptive analysis
    - gc content analysis
  - Protein embedding analysis
  - alignment tools
    - fasta, fastq
    - bowtie2 parameters
    - magicblast
    - bloomine
    - sourmash
    - sam, bam
    - using slurm
      - storage considerations
  - variant calling tools
  - hansel&gretel


- Critical Evaluation (2000)
  - In this chapter (it won’t be chapter 4, but probably chapter 6 or 7 once all the core chapters have been added.) The critical evaluation consists of a discussion, leading to conclusion.  It is an essential part of a master’s degree. It shows that you can not only carry out a substantial piece of work, but that you can reflect on it, and think critically about how you might have done it better. Examiners view the critical evaluation as very important.
  

- Conclusion (500)
  - A brief summary of all that has gone before, including the key results of the project. May also include some directions for future work.
- References
