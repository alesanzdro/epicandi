#!/usr/bin/env python3
import json
import sys

def generate_weights(json_file_path):
    """
    Lee un fichero params.json de funannotate2 train y genera una cadena
    de pesos basada en el rendimiento de cada predictor.
    """
    try:
        with open(json_file_path, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"Error: El fichero '{json_file_path}' no fue encontrado.", file=sys.stderr)
        return None
    except json.JSONDecodeError:
        print(f"Error: El fichero '{json_file_path}' no es un JSON válido.", file=sys.stderr)
        return None

    results = {}
    # Extraer las métricas clave para cada herramienta
    for tool, stats in data.get('abinitio', {}).items():
        train_results = stats.get('train_results', {})
        gene_prec = train_results.get('gene_precision', 0)
        exon_prec = train_results.get('exon_precision', 0)
        results[tool] = {'gene_precision': gene_prec, 'exon_precision': exon_prec}

    # Lógica para asignar pesos
    final_weights = {}
    qualified_predictors = {}

    # 1. Descalificar predictores con muy baja precisión a nivel de gen
    disqualification_threshold = 0.10 # 10%
    for tool, metrics in results.items():
        if metrics['gene_precision'] > disqualification_threshold:
            qualified_predictors[tool] = metrics['exon_precision']
        else:
            final_weights[tool] = 0 # Asignar peso 0 si falla

    # 2. Rankear los predictores calificados por su precisión de exón
    if qualified_predictors:
        # Ordenar de mejor a peor según la precisión del exón
        sorted_tools = sorted(
            qualified_predictors, 
            key=qualified_predictors.get, 
            reverse=True
        )
        
        # 3. Asignar pesos: 2 al mejor, 1 a los demás
        best_tool = sorted_tools[0]
        final_weights[best_tool] = 2
        for tool in sorted_tools[1:]:
            final_weights[tool] = 1

    # 4. Formatear la cadena de salida en un orden consistente
    tool_order = ['augustus', 'genemark', 'glimmerhmm', 'snap']
    weight_parts = []
    for tool in tool_order:
        if tool in final_weights:
            weight_parts.append(f"{tool}:{final_weights[tool]}")
            
    return " ".join(weight_parts)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: ./generate_weights.py <ruta_al_fichero_params.json>", file=sys.stderr)
        sys.exit(1)
        
    json_file = sys.argv[1]
    weights_string = generate_weights(json_file)
    
    if weights_string:
        print(weights_string)
