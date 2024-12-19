### O que fazer
- Utilizar o modelo SIR clássico para ajustar aos dados de infectados nas 12 primeiras semanas de 2021. Os dados estão na planilha I.csv;

- Utilizar o modelo SIR clássico ou o modelo SIRS para ajustar o número de infectados aos dados de 2021 presentes na planilha casos_2021.xlsx. Os dados da planilha devem ser transformados para refletir o número acumulado de casos (infecções). Salvar a planilha como csv para leitura no python.

- Utilizar o modelo SIR clássico ou o SIRS com mortes (conforme modificação feita em sala) para ajustar os dados de infecções e óbitos do ano de 2021 (53 semanas). Os dados de óbitos estão na planilha obtidos_2021.xlsx. Modificar a planilha para armazenar o número acumulado de óbitos ao longo de 2021.

### Considerar, para cada ajuste: 
- tolerância igual a 10-3;
- rand1bin;
- rand2bin;
- best1bin;
- best2bin.

### Para cada ajuste, mostre as seguintes saídas:
- Os valores dos parâmetros obtidos no ajuste e o erro associado.
- Faça um plot mostrando os dados experimentais como pontos e a curva com o resultado do modelo utilizando na simulação os valores obtidos com o ajuste. 