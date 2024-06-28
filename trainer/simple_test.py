import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader

from datasets import ReplicationTimingDataset

class SimpleNN(nn.Module):
    def __init__(self, input_size, hidden_size):
        super(SimpleNN, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(hidden_size, input_size)
    
    def forward(self, x):
        out = self.fc1(x)
        out = self.relu(out)
        out = self.fc2(out)
        return out

input_size = 100  
hidden_size = 50
num_epochs = 5
batch_size = 16
learning_rate = 0.001

train_dataset = ReplicationTimingDataset()
train_loader  = DataLoader(dataset=train_dataset, batch_size=batch_size, shuffle=True)

# 初始化模型、损失函数和优化器
model = SimpleNN(input_size, hidden_size)
criterion = nn.CrossEntropyLoss()
optimizer = optim.Adam(model.parameters(), lr=learning_rate)

# 训练模型
for epoch in range(num_epochs):
    for i, rt in enumerate(train_loader):
        # 将输入和标签转换为变量
        rt = rt.view(-1, input_size)  # 调整输入的形状
        
        # 前向传播
        outputs = model(rt)
        loss = criterion(outputs, rt)
        
        # 反向传播和优化
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        if (i+1) % 100 == 0:
            print(f'Epoch [{epoch+1}/{num_epochs}], Step [{i+1}/{len(train_loader)}], Loss: {loss.item():.4f}')

print('Finished Training')

# 验证代码是否运行成功
# 测试数据集的初始化和加载器
test_dataset = ReplicationTimingDataset()
test_loader = DataLoader(dataset=test_dataset, batch_size=batch_size, shuffle=False)

# 测试模型性能
model.eval()
with torch.no_grad():
    correct = 0
    total = 0
    for rt in test_loader:
        rt = rt.view(-1, input_size)
        outputs = model(rt)
        _, predicted = torch.max(outputs.data, 1)
        total += rt.size(0)
        correct += (predicted == rt).sum().item()
    
    print(f'Accuracy of the model on the test images: {100 * correct / total}%')
