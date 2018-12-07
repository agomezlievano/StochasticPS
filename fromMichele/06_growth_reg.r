library(stargazer);

mydata <- read.csv("eci_growth_regression_table.csv", sep = "\t");

m.2 <- lm(log(growth) ~ log(gdppc), data = mydata);
m.3 <- lm(log(growth) ~ log(gdppc) + eci, data = mydata);
m.4 <- lm(log(growth) ~ log(gdppc) + fitness, data = mydata);
m.5 <- lm(log(growth) ~ log(gdppc) + redfc.fix, data = mydata);
m.6 <- lm(log(growth) ~ log(gdppc) + redfc.var, data = mydata);
m.7 <- lm(log(growth) ~ log(gdppc) + eci + fitness + redfc.fix, data = mydata);
m.8 <- lm(log(growth) ~ log(gdppc) + eci + fitness + redfc.var, data = mydata);

stargazer(m.2, m.3, m.4, m.5, m.6, m.7, m.8, title = "Growth Regressions (Meip).", align = TRUE, df = FALSE, model.numbers = FALSE);

mydata <- read.csv("eci_growth_regression_table_cp.csv", sep = "\t");

m.2 <- lm(growth ~ log(gdppc) + log(nr + 1), data = mydata);
m.3 <- lm(growth ~ log(gdppc) + log(nr + 1) + eci, data = mydata);
m.4 <- lm(growth ~ log(gdppc) + log(nr + 1) + fitness, data = mydata);
m.5 <- lm(growth ~ log(gdppc) + log(nr + 1) + log(redfc.fix), data = mydata);
m.6 <- lm(growth ~ log(gdppc) + log(nr + 1) + log(redfc.var), data = mydata);
m.7 <- lm(growth ~ log(gdppc) + log(nr + 1) + eci + fitness + log(redfc.fix), data = mydata);
m.8 <- lm(growth ~ log(gdppc) + log(nr + 1) + eci + fitness + log(redfc.var), data = mydata);

stargazer(m.2, m.3, m.4, m.5, m.6, m.7, m.8, title = "Growth Regressions (Mcp).", align = TRUE, df = FALSE, model.numbers = FALSE);
