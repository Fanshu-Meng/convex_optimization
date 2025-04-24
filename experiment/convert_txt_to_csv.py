import csv

def txt_to_csv(input_txt_path, output_csv_path):
    """
    将特定格式的txt文档转换为csv文件。
    
    参数:
    - input_txt_path: 输入txt文件路径
    - output_csv_path: 输出csv文件路径
    """
    with open(input_txt_path, 'r', encoding='utf-8') as txt_file:
        lines = txt_file.readlines()
    
    # 定义csv的列名
    headers = [
        "Algorithm", "Test Function", "Optimal Solution", 
        "Function Value", "Loss (L2 Norm)", "Iterations", "Run Time"
    ]
    
    # 存储解析后的数据
    data = []
    current_entry = {}
    
    for line in lines:
        line = line.strip()
        if not line:  # 空行表示一个记录结束
            if current_entry:  # 如果当前记录非空，保存到数据列表
                data.append(current_entry)
                current_entry = {}
        else:
            # 解析每一行
            if line.startswith("Algorithm:"):
                current_entry["Algorithm"] = line.split(":", 1)[1].strip()
            elif line.startswith("Test Function:"):
                current_entry["Test Function"] = line.split(":", 1)[1].strip()
            elif line.startswith("Optimal Solution:"):
                current_entry["Optimal Solution"] = line.split(":", 1)[1].strip()
            elif line.startswith("Function Value:"):
                current_entry["Function Value"] = line.split(":", 1)[1].strip()
            elif line.startswith("Error (L2 Norm):"):
                current_entry["Loss (L2 Norm)"] = line.split(":", 1)[1].strip()
            elif line.startswith("Iterations:"):
                current_entry["Iterations"] = line.split(":", 1)[1].strip()
            elif line.startswith("Run Time:"):
                current_entry["Run Time"] = line.split(":", 1)[1].strip()
    
    # 如果最后一条记录没有被添加
    if current_entry:
        data.append(current_entry)
    
    # 写入csv文件
    with open(output_csv_path, 'w', newline='', encoding='utf-8') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=headers)
        writer.writeheader()
        writer.writerows(data)

txt_to_csv('figures/log.txt', 'figures/log.csv')