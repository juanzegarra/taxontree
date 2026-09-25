

library(dplyr)
library(ggplot2)
library(tydr)



grouped_benchmark <- benchmark_results %>% group_by(program, fasta_file, cores) %>% 
  summarise(mean_c = mean(time_seconds, na.rm = TRUE), 
            std_c = sd(time_seconds, na.rm = TRUE)) %>% 
  ungroup()

grouped_benchmark <- inner_join(grouped_benchmark, len_fasta, by = c("fasta_file" = "fasta"))

library(ggplot2)

program_list <- unique(benchmark_results$program)

for (prog in program_list) {
  
  # 1. Filtra os dados para o programa atual
  df_sub <- benchmark_results %>% filter(program == prog)
  
  # 2. Cria o gráfico com agrupamento por 'cores'
  p <- ggplot(df_sub, aes(
    x = factor(length), 
    y = time_seconds, 
    fill = factor(cores)
  )) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    # Adiciona os pontos das réplicas alinhados com cada caixa respectiva
    geom_point(
      position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
      alpha = 0.6,
      size = 1.5
    ) +
    theme_minimal() +
    labs(
      title = paste("Desempenho por Núcleos (Cores) -", prog),
      x = "Tamanho da Sequência (resíduos)",
      y = "Tempo (segundos)",
      fill = "Nº de Cores"
    ) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold")
    )
  
  # 3. Exibe o gráfico na tela
  print(p)
  
  # 4. Salva a imagem em alta resolução
  ggsave(
    filename = paste0("boxplot_", prog, "_cores.png"),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# extrair maximo e minimo de uso de memoria 
# 3 min e 40s para blast na bioinfo taxontree 0 

