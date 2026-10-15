import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import distance_array
import collections
import pandas as pd
import os
import sys
import time


# 进度条显示函数
def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=50, fill='█'):
    """
    显示进度条

    参数:
        iteration: 当前迭代次数
        total: 总迭代次数
        prefix: 进度条前缀字符串
        suffix: 进度条后缀字符串
        decimals: 进度百分比的小数位数
        length: 进度条长度（字符数）
        fill: 进度条填充字符
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    sys.stdout.write(f'\r{prefix} |{bar}| {percent}% {suffix}')
    sys.stdout.flush()

    # 当完成时打印新行
    if iteration == total:
        print()


# 读取xyz文件
def load_xyz(filename):
    # 使用MDAnalysis加载.xyz文件
    u = mda.Universe(filename, format='XYZ')
    if u.dimensions is None:
        with open(filename, 'r') as file:
            file.readline()
            # 第二行：盒子大小（Lattice信息）
            lattice_line = file.readline().split('"')[1].split()
            lattice = np.array([float(lattice_line[0]), float(lattice_line[4]), float(lattice_line[8]), 90, 90, 90])
        u.dimensions = lattice
    return u


# 修改 calculate_cn 函数以应用 PBC
def calculate_cn(u, cation_types, bond_limits, boundary_condition):
    cations = u.select_atoms("name " + " or name ".join(cation_types))
    oxygens = u.select_atoms("name O")

    cn_stats = {cation: {f"CN{i}": 0 for i in range(16)} for cation in cation_types}
    bonding_cations = {c.index: None for c in cations}
    bonding_oxy = {ox.index: [] for ox in oxygens}

    # 获取模拟盒子尺寸以应用 PBC
    if boundary_condition[0].upper() == 'Y':
        box = u.dimensions
    else:
        box = None

    # 存储每个阳离子的配位数
    cation_cn_dict = {}

    for cation in cations:
        cation_type = cation.name
        bond_limit = bond_limits[cation_type]

        # 获取与该阳离子成键的氧离子并考虑 PBC
        distances = distance_array(oxygens.positions, cation.position, box=box).flatten()
        bonded_oxygens = np.where(distances <= bond_limit)[0]
        cn = len(bonded_oxygens)

        # 记录阳离子的配位数
        cation_cn_dict[cation.index] = cn

        for ox_local_index in bonded_oxygens:
            ox_global_index = oxygens[ox_local_index].index
            bonding_oxy[ox_global_index].append(cation)

        if cn < 16:
            cn_stats[cation_type][f"CN{cn}"] += 1

        bonding_cations[cation.index] = oxygens[bonded_oxygens]

    return cn_stats, bonding_cations, bonding_oxy, cation_cn_dict


# 氧离子种类统计（修复版）
def classify_oxygen(u, bonding_cations, bonding_oxy, forming_B):
    oxygen_types = {
        'Bridge Oxygen': 0,
        'Non Bridge Oxygen': 0,
        'Free Oxygen': 0,
        'Tri Oxygen': 0
    }

    # 用于统计非桥氧中修饰阳离子的元素分布
    non_bridge_modifiers = {}

    # 用于记录已经统计过的修饰阳离子
    counted_modifiers = set()

    oxy_type = dict()  # 记录每个氧离子所属种类
    epoxy = dict()  # 记录与桥氧连接的形成体阳离子
    three_bonding_oxy = dict()  # 记录与三簇氧连接的形成体阳离子

    # 获取所有修饰阳离子
    all_modifiers = u.select_atoms("not name O and not name Si and not name Al and not name B and not name Zr")
    modifier_indices = {atom.index for atom in all_modifiers}

    for ox_index, bonds in bonding_oxy.items():
        forming_cations = 0
        bonding_forming_cations = list()
        cation_Al5or6 = list()

        for cation in bonds:
            # 确保我们处理的是原子对象
            if not hasattr(cation, 'index'):
                continue

            ligancy = len(bonding_cations[cation.index]) if cation.index in bonding_cations else 0

            if cation.name == 'Si':
                forming_cations += 1
                bonding_forming_cations.append(cation)
            elif cation.name == 'Al' and ligancy == 4:
                forming_cations += 1
                bonding_forming_cations.append(cation)
            elif cation.name == 'Zr':
                forming_cations += 1
                bonding_forming_cations.append(cation)
            elif cation.name == 'B' and ligancy in forming_B:
                forming_cations += 1
                bonding_forming_cations.append(cation)
            elif cation.name == 'Al' and ligancy in [5, 6]:
                cation_Al5or6.append(cation)

        if len(cation_Al5or6) + forming_cations == 3:  # 三簇氧
            oxygen_types['Tri Oxygen'] += 1
            oxy_type[ox_index] = "Tri Oxygen"
            three_bonding_oxy[ox_index] = cation_Al5or6 + bonding_forming_cations
        elif forming_cations == 2:  # 桥氧
            oxygen_types['Bridge Oxygen'] += 1
            oxy_type[ox_index] = "Bridge Oxygen"
            epoxy[ox_index] = bonding_forming_cations
        elif forming_cations == 1:  # 非桥氧
            oxygen_types['Non Bridge Oxygen'] += 1
            oxy_type[ox_index] = "Non Bridge Oxygen"

            # 统计非桥氧中的修饰阳离子元素分布
            for cation in bonds:
                # 确保我们处理的是原子对象
                if not hasattr(cation, 'index'):
                    continue

                # 只统计修饰阳离子（非形成体）
                if cation.index in modifier_indices and cation.index not in counted_modifiers:
                    element = cation.name
                    if element in non_bridge_modifiers:
                        non_bridge_modifiers[element] += 1
                    else:
                        non_bridge_modifiers[element] = 1
                    # 标记此修饰阳离子已统计
                    counted_modifiers.add(cation.index)
        else:  # 自由氧
            oxygen_types['Free Oxygen'] += 1
            oxy_type[ox_index] = "Free Oxygen"

    # 修复：使用正确的变量名 oxygen_types
    oxygen_types['Total Oxygen'] = sum(oxygen_types.values())

    # 修复：返回正确的变量名 oxygen_types
    return oxygen_types, oxy_type, epoxy, three_bonding_oxy, non_bridge_modifiers


# Qn统计
def calculate_qn(u, bonding_cations, oxy_type):
    atoms = u.atoms
    qn_stats = dict()

    for cation_index, oxygens in bonding_cations.items():
        bridge_oxygens = 0
        # 遍历该阳离子与之成键的氧离子
        for oxygen in oxygens:
            # 获取该氧离子的成键阳离子数
            if oxy_type[oxygen.index] == "Bridge Oxygen":  # 是桥氧
                bridge_oxygens += 1
        # 根据桥氧的数量对阳离子进行Qn分类
        cation_name = atoms[cation_index].name + str(len(oxygens))  # 阳离子配位数
        if not (cation_name in qn_stats):
            qn_stats[cation_name] = {f"Q{i}": 0 for i in range(6)}

        if bridge_oxygens <= 5:
            qn_stats[cation_name][f"Q{bridge_oxygens}"] += 1

    return collections.OrderedDict(sorted(qn_stats.items()))


# 新增：详细的Qn统计，区分Si和Al连接的Al/Si数量
def calculate_detailed_qn(u, bonding_cations, oxy_type, epoxy):
    """
    详细的Qn统计，区分Si和Al四面体单元连接的Al/Si数量

    参数:
        u: MDAnalysis Universe对象
        bonding_cations: 成键阳离子字典
        oxy_type: 氧类型字典
        epoxy: 桥氧连接的阳离子对

    返回:
        si_detailed_stats: Si的详细统计字典
        al_detailed_stats: Al的详细统计字典
    """
    atoms = u.atoms

    # 初始化统计字典
    si_detailed_stats = {}
    al_detailed_stats = {}

    # 获取所有Si和Al原子
    si_atoms = u.select_atoms("name Si")
    al_atoms = u.select_atoms("name Al")

    # 为每个Si原子统计详细的Qn信息
    for si in si_atoms:
        si_index = si.index
        # 只考虑配位数为4的Si（四面体）
        if len(bonding_cations[si_index]) != 4:
            continue

        # 统计该Si的桥氧数量和连接的Al数量
        bridge_oxy_count = 0
        connected_al_count = 0

        # 遍历该Si连接的氧原子
        for oxygen in bonding_cations[si_index]:
            if oxy_type[oxygen.index] == "Bridge Oxygen":
                bridge_oxy_count += 1
                # 查找该桥氧连接的另一个阳离子
                if oxygen.index in epoxy:
                    for cation in epoxy[oxygen.index]:
                        if cation.index != si_index and cation.name == 'Al' and len(bonding_cations[cation.index]) == 4:
                            connected_al_count += 1

        # 根据桥氧数和连接的Al数构建构型名称
        config_name = f"Q{bridge_oxy_count}({connected_al_count}Al)"

        # 更新统计
        if config_name in si_detailed_stats:
            si_detailed_stats[config_name] += 1
        else:
            si_detailed_stats[config_name] = 1

    # 为每个Al原子统计详细的Qn信息
    for al in al_atoms:
        al_index = al.index
        # 只考虑配位数为4的Al（四面体）
        if len(bonding_cations[al_index]) != 4:
            continue

        # 统计该Al的桥氧数量和连接的Si数量
        bridge_oxy_count = 0
        connected_si_count = 0

        # 遍历该Al连接的氧原子
        for oxygen in bonding_cations[al_index]:
            if oxy_type[oxygen.index] == "Bridge Oxygen":
                bridge_oxy_count += 1
                # 查找该桥氧连接的另一个阳离子
                if oxygen.index in epoxy:
                    for cation in epoxy[oxygen.index]:
                        if cation.index != al_index and cation.name == 'Si' and len(bonding_cations[cation.index]) == 4:
                            connected_si_count += 1

        # 根据桥氧数和连接的Si数构建构型名称
        config_name = f"Q{bridge_oxy_count}({connected_si_count}Si)"

        # 更新统计
        if config_name in al_detailed_stats:
            al_detailed_stats[config_name] += 1
        else:
            al_detailed_stats[config_name] = 1

    # 确保所有可能的构型都存在（即使数量为0）
    possible_si_configs = [f"Q{i}({j}Al)" for i in range(5) for j in range(i + 1)]
    possible_al_configs = [f"Q{i}({j}Si)" for i in range(5) for j in range(i + 1)]

    for config in possible_si_configs:
        if config not in si_detailed_stats:
            si_detailed_stats[config] = 0

    for config in possible_al_configs:
        if config not in al_detailed_stats:
            al_detailed_stats[config] = 0

    # 按构型名称排序
    si_detailed_stats = collections.OrderedDict(sorted(si_detailed_stats.items()))
    al_detailed_stats = collections.OrderedDict(sorted(al_detailed_stats.items()))

    return si_detailed_stats, al_detailed_stats


# NetLink统计
def calculate_net(u, epoxy, bonding_cations):
    atoms = u.atoms
    net_link = dict()

    for _, cations in epoxy.items():
        if len(cations) == 2:
            cation1, cation2 = cations
            cation1_name = cation1.name + str(len(bonding_cations[cation1.index]))
            cation2_name = cation2.name + str(len(bonding_cations[cation2.index]))
            sorted_name = sorted([cation1_name, cation2_name])
            let_name = sorted_name[0] + '-O-' + sorted_name[1]
            if not (let_name in net_link):
                net_link[let_name] = 0
            net_link[let_name] += 1

    return collections.OrderedDict(sorted(net_link.items()))


# 三簇氧统计
def calculate_three_bonding_oxy(u, three_bonding_oxy, bonding_cations):
    atoms = u.atoms
    tbo_stats = dict()

    for _, cations in three_bonding_oxy.items():
        if len(cations) == 3:
            cation1, cation2, cation3 = cations
            cation1_name = cation1.name + str(len(bonding_cations[cation1.index]))
            cation2_name = cation2.name + str(len(bonding_cations[cation2.index]))
            cation3_name = cation3.name + str(len(bonding_cations[cation3.index]))
            sorted_name = sorted([cation1_name, cation2_name, cation3_name])
            tbo_name = sorted_name[0] + '-O-' + sorted_name[1] + '-O-' + sorted_name[2]
            if not (tbo_name in tbo_stats):
                tbo_stats[tbo_name] = 0
            tbo_stats[tbo_name] += 1

    return collections.OrderedDict(sorted(tbo_stats.items()))


# 修改后的分布计算函数：同时计算AlO4、AlO5和三簇氧的分布
def calculate_distributions(u, cation_cn_dict, oxy_type, modifier_type,
                            cutoff_radius=5.0, bin_size=0.02):
    """
    计算以修饰阳离子为中心：
    1. AlO4和AlO5在不同距离处的分布概率
    2. 三簇氧在不同距离处的分布概率

    参数:
    u: MDAnalysis Universe对象
    cation_cn_dict: 阳离子配位数字典 {index: 配位数}
    oxy_type: 氧离子类型字典 {index: 类型}
    modifier_type: 修饰阳离子类型 (如 'Mg')
    cutoff_radius: 截止半径 (Å)
    bin_size: 距离区间大小 (Å)

    返回:
    bin_centers: 距离区间中点列表
    alo4_probabilities: AlO4在各距离区间的概率
    alo5_probabilities: AlO5在各距离区间的概率
    three_oxy_probabilities: 三簇氧在各距离区间的概率
    """
    # 选择所有指定类型的修饰阳离子
    modifiers = u.select_atoms(f"name {modifier_type}")
    if len(modifiers) == 0:
        return np.array([]), np.array([]), np.array([]), np.array([])

    # 选择所有Al原子和所有O原子
    al_atoms = u.select_atoms("name Al")
    all_oxygens = u.select_atoms("name O")

    # 初始化距离区间
    num_bins = int(cutoff_radius / bin_size) + 1
    distance_bins = np.linspace(0, cutoff_radius, num_bins)
    bin_centers = (distance_bins[:-1] + distance_bins[1:]) / 2

    # 初始化计数器
    alo4_counts = np.zeros(len(bin_centers))
    alo5_counts = np.zeros(len(bin_centers))
    three_oxy_counts = np.zeros(len(bin_centers))
    modifier_counts_al = np.zeros(len(bin_centers))  # 用于铝统计的修饰阳离子计数
    modifier_counts_oxy = np.zeros(len(bin_centers))  # 用于氧统计的修饰阳离子计数

    # 获取盒子尺寸以应用PBC
    box = u.dimensions if u.dimensions is not None else None

    # 对每个修饰阳离子计算距离分布
    for modifier in modifiers:
        # 计算修饰阳离子到所有Al原子的距离
        dists_al = distance_array(al_atoms.positions, modifier.position.reshape(1, 3), box=box).flatten()

        # 计算修饰阳离子到所有氧原子的距离
        dists_oxy = distance_array(all_oxygens.positions, modifier.position.reshape(1, 3), box=box).flatten()

        # 对Al原子：统计AlO4和AlO5
        for i, dist in enumerate(dists_al):
            if dist <= cutoff_radius:
                # 找到对应的距离区间
                bin_idx = int(dist / bin_size)
                if bin_idx < len(bin_centers):
                    # 获取Al原子的配位数
                    al_idx = al_atoms[i].index
                    cn = cation_cn_dict.get(al_idx, 0)

                    # 根据配位数分类计数
                    if cn == 4:
                        alo4_counts[bin_idx] += 1
                    elif cn == 5:
                        alo5_counts[bin_idx] += 1

                    # 记录修饰阳离子的数量（用于归一化）
                    modifier_counts_al[bin_idx] += 1

        # 对氧原子：统计三簇氧
        for i, dist in enumerate(dists_oxy):
            if dist <= cutoff_radius:
                # 找到对应的距离区间
                bin_idx = int(dist / bin_size)
                if bin_idx < len(bin_centers):
                    ox_idx = all_oxygens[i].index
                    # 检查是否为三簇氧
                    if oxy_type.get(ox_idx) == "Tri Oxygen":
                        three_oxy_counts[bin_idx] += 1

                    # 记录修饰阳离子的数量（用于归一化）
                    modifier_counts_oxy[bin_idx] += 1

    # 计算概率（归一化）
    alo4_probabilities = np.zeros(len(bin_centers))
    alo5_probabilities = np.zeros(len(bin_centers))
    three_oxy_probabilities = np.zeros(len(bin_centers))

    for i in range(len(bin_centers)):
        if modifier_counts_al[i] > 0:
            alo4_probabilities[i] = alo4_counts[i] / modifier_counts_al[i]
            alo5_probabilities[i] = alo5_counts[i] / modifier_counts_al[i]

        if modifier_counts_oxy[i] > 0:
            three_oxy_probabilities[i] = three_oxy_counts[i] / modifier_counts_oxy[i]

    return bin_centers, alo4_probabilities, alo5_probabilities, three_oxy_probabilities


# 修改功能：细化AlO4/5/6中三簇氧的数量统计
def calculate_tri_oxygen_in_alo(u, bonding_cations, cation_cn_dict, oxy_type):
    """
    细化统计AlO4、AlO5和AlO6中三簇氧的连接数量

    参数:
    u: MDAnalysis Universe对象
    bonding_cations: 成键信息字典 {阳离子index: 连接的氧原子列表}
    cation_cn_dict: 阳离子配位数字典 {index: 配位数}
    oxy_type: 氧离子类型字典 {index: 类型}

    返回:
    包含详细统计结果的字典
    """
    # 选择所有铝原子
    al_atoms = u.select_atoms("name Al")

    # 初始化统计字典
    # 键是配位数类型，值是一个字典，键是连接的三簇氧数量，值是铝原子数量
    tri_oxygen_stats = {
        'AlO4': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0},
        'AlO5': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0},
        'AlO6': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0}
    }

    # 用于计算平均值的辅助字典
    avg_stats = {
        'AlO4': {'total_tri_oxy': 0, 'count': 0, 'avg_tri_per_alo': 0.0},
        'AlO5': {'total_tri_oxy': 0, 'count': 0, 'avg_tri_per_alo': 0.0},
        'AlO6': {'total_tri_oxy': 0, 'count': 0, 'avg_tri_per_alo': 0.0}
    }

    # 遍历所有铝原子
    for al in al_atoms:
        al_index = al.index
        cn = cation_cn_dict.get(al_index, 0)

        # 只处理配位数为4、5、6的铝原子
        if cn in [4, 5, 6]:
            # 获取该铝原子连接的氧原子
            bonded_oxy = bonding_cations.get(al_index, [])
            tri_oxy_count = 0

            # 统计三簇氧的数量
            for ox in bonded_oxy:
                if oxy_type.get(ox.index) == "Tri Oxygen":
                    tri_oxy_count += 1

            # 根据配位数分类统计
            if cn == 4:
                # 确保不会超出范围
                tri_oxy_count = min(tri_oxy_count, 4)
                tri_oxygen_stats['AlO4'][tri_oxy_count] += 1
                avg_stats['AlO4']['total_tri_oxy'] += tri_oxy_count
                avg_stats['AlO4']['count'] += 1
            elif cn == 5:
                tri_oxy_count = min(tri_oxy_count, 5)
                tri_oxygen_stats['AlO5'][tri_oxy_count] += 1
                avg_stats['AlO5']['total_tri_oxy'] += tri_oxy_count
                avg_stats['AlO5']['count'] += 1
            elif cn == 6:
                tri_oxy_count = min(tri_oxy_count, 6)
                tri_oxygen_stats['AlO6'][tri_oxy_count] += 1
                avg_stats['AlO6']['total_tri_oxy'] += tri_oxy_count
                avg_stats['AlO6']['count'] += 1

    # 计算平均每个AlO的三簇氧数量
    for key in avg_stats:
        if avg_stats[key]['count'] > 0:
            avg_stats[key]['avg_tri_per_alo'] = (
                    avg_stats[key]['total_tri_oxy'] / avg_stats[key]['count']
            )

    # 返回详细统计和平均值统计
    return tri_oxygen_stats, avg_stats


# 修改保存函数以包含详细的三簇氧连接统计和详细的Qn统计
def save_to_excel(cn_stats, oxygen_stats, qn_stats, net_link, tbo_stats,
                  non_bridge_modifiers, dist_data, tri_oxygen_global, avg_tri_global, total_files,
                  si_detailed_global, al_detailed_global, output_filename):
    # 1 配位数统计
    df_cn = pd.DataFrame(cn_stats)
    df_cn.loc['Total'] = df_cn.sum(axis=0)
    df_cn['Total'] = df_cn.sum(axis=1)
    df_cn = df_cn.T

    # 2 氧类型统计
    df_oxygen = pd.DataFrame(oxygen_stats, index=['number']).T

    # 3 Qn统计
    df_qn = pd.DataFrame(qn_stats)
    df_qn.loc['Total'] = df_qn.sum(axis=0)
    df_qn['Total'] = df_qn.sum(axis=1)
    df_qn = df_qn.T

    # 4 桥氧网络统计
    df_net_link = pd.DataFrame(net_link, index=['number'])
    df_net_link['Total'] = df_net_link.sum(axis=1)
    df_net_link = df_net_link.T

    # 5 三簇氧统计
    df_tbo = pd.DataFrame(tbo_stats, index=['number'])
    df_tbo['Total'] = df_tbo.sum(axis=1)
    df_tbo = df_tbo.T

    # 6 非桥氧修饰阳离子统计
    df_modifiers = pd.DataFrame(non_bridge_modifiers, index=['count']).T
    df_modifiers.loc['Total'] = df_modifiers.sum()

    # 7 分布统计（AlO4, AlO5, 三簇氧）
    dist_dfs = {}

    # 8 AlO4/5/6中三簇氧详细统计
    # 创建一个空的DataFrame用于存储详细统计
    detailed_stats = []

    for alo_type in ['AlO4', 'AlO5', 'AlO6']:
        # 计算每个文件的平均值
        avg_count = avg_tri_global[alo_type]['count'] / total_files
        avg_total_tri_oxy = avg_tri_global[alo_type]['total_tri_oxy'] / total_files

        # 添加平均值行
        detailed_stats.append({
            'AlO Type': alo_type,
            'Tri Oxygen Count': 'Average per File',
            'AlO Count': avg_count,
            'Tri Oxygen Total': avg_total_tri_oxy,
            'Avg Tri per AlO': avg_tri_global[alo_type]['avg_tri_per_alo']
        })

        # 添加每种三簇氧连接数量的统计
        max_tri = 6 if alo_type == 'AlO6' else (5 if alo_type == 'AlO5' else 4)
        for tri_count in range(max_tri + 1):
            count = tri_oxygen_global[alo_type][tri_count] / total_files
            detailed_stats.append({
                'AlO Type': alo_type,
                'Tri Oxygen Count': f'{tri_count}',
                'AlO Count': count,
                'Tri Oxygen Total': count * tri_count,
                'Avg Tri per AlO': tri_count
            })

    df_detailed = pd.DataFrame(detailed_stats)

    # 9 详细的Qn统计
    # 创建详细的Qn统计DataFrame
    detailed_qn_data = []

    # 获取所有可能的Si和Al构型
    si_configs = sorted(si_detailed_global.keys())
    al_configs = sorted(al_detailed_global.keys())

    # 计算Si和Al的总数
    si_total = sum(si_detailed_global.values()) / total_files
    al_total = sum(al_detailed_global.values()) / total_files

    # 确保两个列表长度相同，用空字符串填充
    max_len = max(len(si_configs), len(al_configs))

    for i in range(max_len):
        si_config = si_configs[i] if i < len(si_configs) else ""
        si_count = si_detailed_global[si_config] / total_files if i < len(si_configs) else ""
        si_percentage = (si_count / si_total * 100) if si_count != "" and si_total > 0 else ""

        al_config = al_configs[i] if i < len(al_configs) else ""
        al_count = al_detailed_global[al_config] / total_files if i < len(al_configs) else ""
        al_percentage = (al_count / al_total * 100) if al_count != "" and al_total > 0 else ""

        detailed_qn_data.append({
            'Si Configuration': si_config,
            'Si Count': si_count,
            'Si Percentage (%)': si_percentage,
            '': '',  # 空列作为分隔
            'Al Configuration': al_config,
            'Al Count': al_count,
            'Al Percentage (%)': al_percentage
        })

    # 添加总数行
    detailed_qn_data.append({
        'Si Configuration': 'Total',
        'Si Count': si_total,
        'Si Percentage (%)': 100.0 if si_total > 0 else 0.0,
        '': '',
        'Al Configuration': 'Total',
        'Al Count': al_total,
        'Al Percentage (%)': 100.0 if al_total > 0 else 0.0
    })

    df_detailed_qn = pd.DataFrame(detailed_qn_data)

    with pd.ExcelWriter(output_filename) as writer:
        df_cn.to_excel(writer, sheet_name="Coordination Numbers")
        df_oxygen.to_excel(writer, sheet_name="Oxygen Types")
        df_qn.to_excel(writer, sheet_name="Qn Statistics")
        df_detailed_qn.to_excel(writer, sheet_name="Detailed Qn Statistics", index=False)
        df_net_link.to_excel(writer, sheet_name="Net Link")
        df_tbo.to_excel(writer, sheet_name="Three Bonding Oxygen Statistics")
        df_modifiers.to_excel(writer, sheet_name="Non-Bridge Modifiers")
        df_detailed.to_excel(writer, sheet_name="Detailed Tri Oxygen in AlO", index=False)

        # 为每种修饰体原子创建单独的工作表
        for modifier_type, data in dist_data.items():
            if data['distance'] is not None and len(data['distance']) > 0:
                # 确保所有数组长度一致
                min_len = min(len(data['distance']),
                              len(data['alo4_prob']),
                              len(data['alo5_prob']),
                              len(data['three_oxy_prob']))

                # 创建数据框
                df_dist = pd.DataFrame({
                    'Distance (Å)': data['distance'][:min_len],
                    'AlO4 Probability': data['alo4_prob'][:min_len],
                    'AlO5 Probability': data['alo5_prob'][:min_len],
                    'Three Oxygen Probability': data['three_oxy_prob'][:min_len]
                })

                # 保存到单独的工作表
                sheet_name = f"Distributions ({modifier_type})"
                df_dist.to_excel(writer, sheet_name=sheet_name, index=False)


def main():
    # 用户输入
    folder_path = input("输入.xyz文件所在文件夹的路径: ")
    cation_types = input("输入需要统计的阳离子种类 (写用空格分隔): ").split()
    limits = input("对应阳离子的成键上限 (写用空格分隔): ").split()
    if len(cation_types) != len(limits):
        exit("阳离子种类数目和成键上限数目输入数据不匹配")
    forming_B = [4]
    if ('B' in cation_types or 'b' in cation_types) and input("三配位的B是否作为形成体: (yes/no)")[0].lower() == 'y':
        forming_B.append(3)
    boundary_periodicity = input("周期性边界条件: (Yes/No) ")
    bond_limits = {cation_types[i]: float(limits[i]) for i in range(len(limits))}

    # 获取多种修饰阳离子类型
    modifier_types = input("请输入要分析的修饰阳离子类型 (用空格分隔，如 Mg Ca): ").split()
    if not modifier_types:
        print("未指定修饰阳离子类型，跳过分布统计")

    # 全局统计变量
    cn_global = dict()
    oxygen_global = dict()
    qn_global = dict()
    net_global = dict()
    tbo_global = dict()

    # 全局非桥氧修饰阳离子统计
    non_bridge_modifiers_global = {}

    # 全局分布统计，每种修饰体原子单独存储
    dist_global = {}
    for modifier_type in modifier_types:
        dist_global[modifier_type] = {
            'distance': None,
            'alo4_prob': np.array([]),
            'alo5_prob': np.array([]),
            'three_oxy_prob': np.array([]),
            'count': 0
        }

    # 全局AlO中三簇氧详细统计
    tri_oxygen_global = {
        'AlO4': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0},
        'AlO5': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0},
        'AlO6': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0}
    }

    # 全局AlO中三簇氧平均值统计
    avg_tri_global = {
        'AlO4': {'total_tri_oxy': 0, 'count': 0, 'avg_tri_per_alo': 0.0},
        'AlO5': {'total_tri_oxy': 0, 'count': 0, 'avg_tri_per_alo': 0.0},
        'AlO6': {'total_tri_oxy': 0, 'count': 0, 'avg_tri_per_alo': 0.0}
    }

    # 全局详细的Qn统计
    si_detailed_global = {}
    al_detailed_global = {}

    # 获取所有.xyz文件列表
    xyz_files = [f for f in os.listdir(folder_path) if f.endswith(".xyz")]
    total_files = len(xyz_files)

    if total_files == 0:
        print("没有找到.xyz文件！")
        return

    print(f"找到 {total_files} 个.xyz文件")
    print(f"将分析 {len(modifier_types)} 种修饰阳离子: {', '.join(modifier_types)}")
    print("开始处理文件...")

    start_time = time.time()

    # 遍历文件夹中的所有xyz文件
    for idx, filename in enumerate(xyz_files):
        # 显示进度条
        print_progress_bar(idx, total_files, prefix='进度:', suffix=f'处理文件 {idx + 1}/{total_files}')

        xyz_file = os.path.join(folder_path, filename)
        u = load_xyz(xyz_file)

        # 计算配位数并记录成键信息
        cn_stats, bonding_cations, bonding_oxy, cation_cn_dict = calculate_cn(
            u, cation_types, bond_limits, boundary_periodicity)

        # 累加配位数统计
        for key, value in cn_stats.items():
            if key in cn_global:
                for k, v in value.items():
                    cn_global[key][k] += v
            else:
                cn_global[key] = value

        # 氧离子分类统计
        oxygen_stats, oxy_type, epoxy, three_bonding_oxy, non_bridge_modifiers = classify_oxygen(
            u, bonding_cations, bonding_oxy, forming_B)

        # 累加氧类型统计
        for key, value in oxygen_stats.items():
            if key in oxygen_global:
                oxygen_global[key] += value
            else:
                oxygen_global[key] = value

        # 累加非桥氧修饰阳离子统计
        for element, count in non_bridge_modifiers.items():
            if element in non_bridge_modifiers_global:
                non_bridge_modifiers_global[element] += count
            else:
                non_bridge_modifiers_global[element] = count

        # Qn统计
        qn_stats = calculate_qn(u, bonding_cations, oxy_type)

        # 累加Qn统计
        for key, value in qn_stats.items():
            if key in qn_global:
                for k, v in value.items():
                    qn_global[key][k] += v
            else:
                qn_global[key] = value

        # 详细的Qn统计
        si_detailed_stats, al_detailed_stats = calculate_detailed_qn(u, bonding_cations, oxy_type, epoxy)

        # 累加详细的Qn统计
        for config, count in si_detailed_stats.items():
            if config in si_detailed_global:
                si_detailed_global[config] += count
            else:
                si_detailed_global[config] = count

        for config, count in al_detailed_stats.items():
            if config in al_detailed_global:
                al_detailed_global[config] += count
            else:
                al_detailed_global[config] = count

        # 桥氧网络统计
        net_link = calculate_net(u, epoxy, bonding_cations)

        # 累加桥氧网络统计
        for key, value in net_link.items():
            if key in net_global:
                net_global[key] += value
            else:
                net_global[key] = value

        # 三簇氧统计
        tbo_stats = calculate_three_bonding_oxy(u, three_bonding_oxy, bonding_cations)

        # 累加三簇氧统计
        for key, value in tbo_stats.items():
            if key in tbo_global:
                tbo_global[key] += value
            else:
                tbo_global[key] = value

        # 计算AlO4/5/6中三簇氧的数量（详细统计）
        tri_oxygen_stats, avg_stats = calculate_tri_oxygen_in_alo(u, bonding_cations, cation_cn_dict, oxy_type)

        # 累加AlO中三簇氧详细统计
        for alo_type in ['AlO4', 'AlO5', 'AlO6']:
            for tri_count, count in tri_oxygen_stats[alo_type].items():
                tri_oxygen_global[alo_type][tri_count] += count

            # 累加平均值统计
            avg_tri_global[alo_type]['total_tri_oxy'] += avg_stats[alo_type]['total_tri_oxy']
            avg_tri_global[alo_type]['count'] += avg_stats[alo_type]['count']

        # 对每种修饰体原子分别计算分布（AlO4, AlO5, 三簇氧）
        for modifier_type in modifier_types:
            # 计算当前修饰体原子的分布
            distance_bins, alo4_probs, alo5_probs, three_oxy_probs = calculate_distributions(
                u, cation_cn_dict, oxy_type, modifier_type)

            # 获取当前修饰体原子的统计对象
            dist_data = dist_global[modifier_type]

            # 初始化距离数组（只在第一次运行时）
            if dist_data['distance'] is None and len(distance_bins) > 0:
                dist_data['distance'] = distance_bins

            # 累加分布统计
            if (len(alo4_probs) > 0 and len(alo5_probs) > 0 and len(three_oxy_probs) > 0):
                # 如果是第一次，初始化数组
                if dist_data['alo4_prob'].size == 0:
                    dist_data['alo4_prob'] = np.zeros_like(alo4_probs)
                    dist_data['alo5_prob'] = np.zeros_like(alo5_probs)
                    dist_data['three_oxy_prob'] = np.zeros_like(three_oxy_probs)

                # 确保数组长度一致
                min_len = min(len(alo4_probs), len(dist_data['alo4_prob']),
                              len(alo5_probs), len(dist_data['alo5_prob']),
                              len(three_oxy_probs), len(dist_data['three_oxy_prob']))

                dist_data['alo4_prob'][:min_len] += alo4_probs[:min_len]
                dist_data['alo5_prob'][:min_len] += alo5_probs[:min_len]
                dist_data['three_oxy_prob'][:min_len] += three_oxy_probs[:min_len]
                dist_data['count'] += 1

    # 完成进度条
    print_progress_bar(total_files, total_files, prefix='进度:', suffix=f'处理完成 {total_files}/{total_files}')

    # 计算平均值（每个文件）
    for key, value in cn_global.items():
        for k in value:
            cn_global[key][k] = value[k] / total_files

    for key in oxygen_global:
        oxygen_global[key] = oxygen_global[key] / total_files

    for key, value in qn_global.items():
        for k in value:
            qn_global[key][k] = value[k] / total_files

    for key in net_global:
        net_global[key] = net_global[key] / total_files

    for key in tbo_global:
        tbo_global[key] = tbo_global[key] / total_files

    # 计算非桥氧修饰阳离子平均值
    for element in non_bridge_modifiers_global:
        non_bridge_modifiers_global[element] = non_bridge_modifiers_global[element] / total_files

    # 计算分布平均值
    for modifier_type in modifier_types:
        dist_data = dist_global[modifier_type]
        if dist_data['count'] > 0 and dist_data['alo4_prob'].size > 0:
            dist_data['alo4_prob'] /= dist_data['count']
            dist_data['alo5_prob'] /= dist_data['count']
            dist_data['three_oxy_prob'] /= dist_data['count']

    # 计算AlO中三簇氧的平均值
    for alo_type in ['AlO4', 'AlO5', 'AlO6']:
        if avg_tri_global[alo_type]['count'] > 0:
            avg_tri_global[alo_type]['avg_tri_per_alo'] = (
                    avg_tri_global[alo_type]['total_tri_oxy'] /
                    avg_tri_global[alo_type]['count']
            )
        else:
            avg_tri_global[alo_type]['avg_tri_per_alo'] = 0.0

    # 保存结果
    output_filename = "multifile_cn_and_oxygen_statistics.xlsx"
    save_to_excel(cn_global, oxygen_global, qn_global, net_global, tbo_global,
                  non_bridge_modifiers_global, dist_global, tri_oxygen_global, avg_tri_global, total_files,
                  si_detailed_global, al_detailed_global, output_filename)


    # 计算并显示总运行时间
    end_time = time.time()
    total_time = end_time - start_time
    print(f"处理完成！总共耗时: {total_time:.2f} 秒")
    print(f"结果已保存至 {output_filename}")


if __name__ == "__main__":
    main()