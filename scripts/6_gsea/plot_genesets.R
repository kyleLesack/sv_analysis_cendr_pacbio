library(ggplot2)
set.seed(100)

df2=data.frame(class=sample(c("BP","MF","CC"), 30, replace=T), "GO Term"=paste0("term_", sample(1:30, 30, replace = F)), count=sample(10:30, 30, replace = T), comparison=sample(c("A","B","C"), 30, replace = T))

ggplot(df2, aes(count, GO.Term, fill = class)) +
	geom_bar(stat = "identity") +
	facet_wrap( ~ comparison, scales = "free", nrow = 3) +
	geom_text(
		aes(label = count),
		color = "black",
		hjust = -0.1,
		size = 4,
		position = position_dodge(0.9)
	) +
	theme(
		legend.position = "bottom",
		panel.grid = element_blank(),
		axis.text.x = element_blank(),
		axis.ticks = element_blank(),
		axis.title.y = element_blank(),
		strip.text.x = element_text(size = 14, face = "bold"),
		strip.background = element_blank()
	)

