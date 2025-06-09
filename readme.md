# Prova 1  Tópicos em Computação II
Este repositório contém a implementação de uma estrutura de dados Kd-tree como parte da avaliação da disciplina de Tópicos em Computação II.

#### Tarefa A1
Completar a definição da classe de nós da árvore, KdTree<D, R, A>::Node. Conforme discutido em sala, a classe pode representar tanto ramos quanto folhas da árvore ou, alternativamente, ser derivada em duas classes, uma que representa ramos e outra que representa folhas. Todo nó deve conhecer sua AABB e seu nível na árvore. Um ramo possui um ponteiro para seus dois nós filhos (que podem ser alocados sequencialmente em memória dinâmica). Uma folha contém um índice usado para determinar qual é o primeiro ponto do subconjunto de pontos do nó e um contador do número de pontos.

#### Tarefa A2
Completar a implementação do construtor KdTree(A&& points, const Params& params);. A execução do método deve resultar na árvore refinada de acordo com o critério de subdivisão descrito anteriormente. Acrescente na classe a definição de quaisquer novos membros que julgar necessários.

#### Tarefa A3
Completar a classe KNN e implementar os métodos de KNN e de busca de vizinhos por raio declarados em KdTree.

#### Tarefa A4
Escrever funções de teste (que devem ser chamadas na função principal em Main.cpp) de todas as funcionalidades da classe KdTree. O programa deve funcionar para R²
(para isso, defina classes de pontos e AABBs no plano) e R³ 
e para qualquer tipo de número real. Como tipos de conjuntos de pontos, considere um vetor de pontos simples(como definido em p1/Utils.h) e também o sistema de partículas objeto do Exercício 3. Para o último, dadoque toda partícula tem uma cor como um de seus atributos, escreva uma função de filtro para descartar do KNNe da busca por raio todas as partículas que não tiverem uma determinada cor. Para testar o resultado dasbuscas, escreva uma função que imprima, na saída padrão, o ponto dado e a distância, o índice em A e ascoordenadas de todos os pontos mais próximos. Teste as buscas para valores distintos de k e radius econjuntos com número variável de pontos (incluindo conjuntos com centenas de milhares de pontos).

### Como Compilar e Executar

### Compilação
```
make
```

### Execução
```
./bin/p1kd
```

### Compilação e Execução
```
make run
```