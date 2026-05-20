import arcpy
import os

# Parámetros
red_drenaje = arcpy.GetParameterAsText(0)
flow_dir_raster = arcpy.GetParameterAsText(1)
limite_predio = arcpy.GetParameterAsText(2)
output_fc = arcpy.GetParameterAsText(3)

# Desplazamientos D8
d8 = {
    1: (1, 0), 2: (1, -1), 4: (0, -1), 8: (-1, -1),
    16: (-1, 0), 32: (-1, 1), 64: (0, 1), 128: (1, 1)
}

desc_raster = arcpy.Describe(flow_dir_raster)
cell_size = desc_raster.meanCellWidth
sr = desc_raster.spatialReference

# ═══════════════════════════════════════════════════════
# FASE 1: Identificar puntos de cruce con el borde
# ═══════════════════════════════════════════════════════

arcpy.AddMessage("Fase 1: Identificando puntos de cruce...")

borde = os.path.join("in_memory", "borde_predio")
arcpy.management.PolygonToLine(limite_predio, borde)

puntos_cruce = os.path.join("in_memory", "puntos_cruce")
arcpy.analysis.Intersect([red_drenaje, borde], puntos_cruce, output_type="POINT")

count = int(arcpy.management.GetCount(puntos_cruce)[0])
arcpy.AddMessage(f"   Puntos de cruce: {count}")

if count == 0:
    arcpy.AddError("No se encontraron intersecciones.")
    raise SystemExit

puntos_singles = os.path.join("in_memory", "puntos_singles")
arcpy.management.MultipartToSinglepart(puntos_cruce, puntos_singles)

# ═══════════════════════════════════════════════════════
# FASE 2: Filtrar solo puntos de SALIDA
# ═══════════════════════════════════════════════════════

arcpy.AddMessage("Fase 2: Identificando puntos de salida...")

predio_geom = None
with arcpy.da.SearchCursor(limite_predio, ["SHAPE@"]) as cursor:
    for row in cursor:
        predio_geom = row[0]
        break

puntos_salida = []

with arcpy.da.SearchCursor(puntos_singles, ["SHAPE@XY"]) as cursor:
    for row in cursor:
        xy = row[0]
        if xy is None:
            continue
        x, y = xy[0], xy[1]
        
        try:
            result = arcpy.management.GetCellValue(flow_dir_raster, f"{x} {y}")
            flow_code = int(result.getOutput(0))
        except:
            continue
        
        if flow_code not in d8:
            continue
        
        dx, dy = d8[flow_code]
        downstream_point = arcpy.PointGeometry(
            arcpy.Point(x + dx * cell_size, y + dy * cell_size), sr
        )
        
        if not predio_geom.contains(downstream_point):
            puntos_salida.append((x, y))
            arcpy.AddMessage(f"   Salida: X={x:.1f}, Y={y:.1f}")

arcpy.AddMessage(f"   Total de salidas: {len(puntos_salida)}")

if len(puntos_salida) == 0:
    arcpy.AddError("No se identificaron puntos de salida.")
    raise SystemExit

# ═══════════════════════════════════════════════════════
# FASE 3: Leer topología de la red
# ═══════════════════════════════════════════════════════

arcpy.AddMessage("Fase 3: Leyendo topologia de la red...")

campos = [f.name for f in arcpy.ListFields(red_drenaje)]
from_field = [c for c in campos if c.lower() == "from_node"][0]
to_field = [c for c in campos if c.lower() == "to_node"][0]
arc_field = [c for c in campos if c.lower() == "arcid"][0] if "arcid" in [c.lower() for c in campos] else "OID@"

# segmentos: arcid -> (from_node, to_node, shape)
segmentos = {}
# hijos_de: to_node -> lista de arcids cuyo FROM_NODE es ese nodo (aguas abajo)
# padres_de: to_node -> lista de arcids cuyo TO_NODE es ese nodo (aguas arriba)
aguas_arriba_de = {}  # para un segmento, cuáles son sus padres (segmentos aguas arriba)

# Indexar segmentos por su TO_NODE
segmentos_por_to_node = {}
segmentos_por_from_node = {}

with arcpy.da.SearchCursor(red_drenaje, [arc_field, from_field, to_field, "SHAPE@"]) as cursor:
    for row in cursor:
        arcid, from_node, to_node, shape = row
        segmentos[arcid] = (from_node, to_node, shape)
        segmentos_por_to_node.setdefault(to_node, []).append(arcid)
        segmentos_por_from_node.setdefault(from_node, []).append(arcid)

arcpy.AddMessage(f"   Segmentos de red: {len(segmentos)}")

# ═══════════════════════════════════════════════════════
# FASE 4: Para cada pour point, encontrar segmento y rastrear aguas arriba
# ═══════════════════════════════════════════════════════

arcpy.AddMessage("Fase 4: Rastreando aguas arriba para cada salida...")

tolerancia = cell_size * 3

def obtener_segmentos_aguas_arriba(arcid_inicial):
    """Devuelve el conjunto de todos los segmentos aguas arriba de un arcid dado."""
    visitados = set()
    pila = [arcid_inicial]
    
    while pila:
        actual = pila.pop()
        if actual in visitados:
            continue
        visitados.add(actual)
        
        # Encontrar padres: segmentos cuyo TO_NODE es el FROM_NODE del actual
        from_node_actual = segmentos[actual][0]
        padres = segmentos_por_to_node.get(from_node_actual, [])
        
        for padre in padres:
            if padre != actual and padre not in visitados:
                pila.append(padre)
    
    return visitados

# Para cada pour point, encontrar su segmento y calcular su árbol aguas arriba
pour_data = []  # lista de (x, y, arcid_segmento, conjunto_aguas_arriba)

for x, y in puntos_salida:
    pt_geom = arcpy.PointGeometry(arcpy.Point(x, y), sr)
    
    segmento_mas_cercano = None
    distancia_min = float('inf')
    
    for arcid, (_, _, shape) in segmentos.items():
        dist = pt_geom.distanceTo(shape)
        if dist < distancia_min:
            distancia_min = dist
            segmento_mas_cercano = arcid
    
    if distancia_min > tolerancia:
        arcpy.AddWarning(f"   Punto ({x:.1f}, {y:.1f}) lejos de la red ({distancia_min:.2f}m).")
    
    aguas_arriba = obtener_segmentos_aguas_arriba(segmento_mas_cercano)
    aguas_arriba.add(segmento_mas_cercano)  # Incluir el propio segmento
    
    pour_data.append((x, y, segmento_mas_cercano, aguas_arriba))
    arcpy.AddMessage(f"   Salida ({x:.1f}, {y:.1f}): segmento={segmento_mas_cercano}, aguas_arriba={len(aguas_arriba)} segmentos")

# ═══════════════════════════════════════════════════════
# FASE 5: Agrupar pour points por solapamiento de árboles aguas arriba
# ═══════════════════════════════════════════════════════

arcpy.AddMessage("Fase 5: Agrupando por arroyo (solapamiento aguas arriba)...")

# Dos pour points pertenecen al mismo arroyo si sus árboles aguas arriba se solapan
# (comparten al menos un segmento)

n = len(pour_data)
parent = list(range(n))

def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

def union(x, y):
    rx, ry = find(x), find(y)
    if rx != ry:
        parent[rx] = ry

for i in range(n):
    for j in range(i + 1, n):
        _, _, _, arriba_i = pour_data[i]
        _, _, _, arriba_j = pour_data[j]
        # Si sus árboles aguas arriba comparten al menos un segmento, son el mismo arroyo
        if arriba_i & arriba_j:
            union(i, j)

# Normalizar IDs a 1, 2, 3...
grupos_unicos = {}
next_id = 1
arroyo_ids = []

for i in range(n):
    raiz = find(i)
    if raiz not in grupos_unicos:
        grupos_unicos[raiz] = next_id
        next_id += 1
    arroyo_ids.append(grupos_unicos[raiz])

arcpy.AddMessage(f"   Arroyos identificados: {len(grupos_unicos)}")

# ═══════════════════════════════════════════════════════
# FASE 6: Escribir feature class de salida
# ═══════════════════════════════════════════════════════

arcpy.AddMessage("Fase 6: Escribiendo pour points...")

out_path = os.path.dirname(output_fc)
out_name = os.path.basename(output_fc)

if arcpy.Exists(output_fc):
    arcpy.management.Delete(output_fc)

arcpy.management.CreateFeatureclass(out_path, out_name, "POINT", spatial_reference=sr)
arcpy.management.AddField(output_fc, "PourID", "LONG")
arcpy.management.AddField(output_fc, "ArroyoID", "LONG")

with arcpy.da.InsertCursor(output_fc, ["SHAPE@", "PourID", "ArroyoID"]) as cursor:
    for i, (x, y, _, _) in enumerate(pour_data):
        arroyo_id = arroyo_ids[i]
        pt_geom = arcpy.PointGeometry(arcpy.Point(x, y), sr)
        cursor.insertRow([pt_geom, i + 1, arroyo_id])
        arcpy.AddMessage(f"   Pour Point {i+1}: X={x:.1f}, Y={y:.1f}, Arroyo={arroyo_id}")

# Limpiar
arcpy.management.Delete(borde)
arcpy.management.Delete(puntos_cruce)
arcpy.management.Delete(puntos_singles)

arcpy.AddMessage(f"\nPour points: {len(pour_data)} puntos en {len(grupos_unicos)} arroyos")